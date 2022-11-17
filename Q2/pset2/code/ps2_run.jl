#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 2 (Q2)
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps2_model.jl        creates functions for models
    *    ps2_figures.jl      creates tables
=#


#======================================================#
#       SET UP
#------------------------------------------------------#

# define directory paths
#   NOTE: in visual studio code, make sure you have the computationalF22 folder opened, so that
#   pwd() automatically returns the file path to that folder
if pwd() == "C:\\Users\\79267"
    cd("C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22\\Q2")
end

root     = joinpath(pwd(), "pset2")
codepath = joinpath(root, "code")
tabpath  = joinpath(root, "tables")
# import packages used to run the model

# Pkg.add("StatFiles")
using Parameters, DataFrames, StatFiles, ForwardDiff, Latexify, CSV, Distributions

# import model functions
include("ps2_model.jl")
include("ps2_tables.jl")

#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Prepare the data
# Read the dataset from
if pwd() == "C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22\\Q2"
    filename = "C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22\\Q2\\pset2\\code\\Mortgage_performance_data.dta"
    file_weigths_dim1 = "C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22\\Q2\\pset2\\code\\nodes_weights_dim1.csv"
    file_weigths_dim2 = "C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22\\Q2\\pset2\\code\\nodes_weights_dim2.csv"
end

function prepare_data(; filename::String)
    df = DataFrame(load(filename))

    # Observed characteristics - time-invariant
    X = copy(df)
    X = select!(X,  ["score_0", "rate_spread", "i_large_loan", "i_medium_loan", "i_refinance", "age_r",  "cltv",
                "dti", "cu", "first_mort_r", "i_FHA",
                 "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"])
#    X = Matrix(X)
    # Add the column for constants
    X_add = ones(length(X[:, 1]))
    X = hcat(X, X_add)

    # Observed characteristics - time-invariant
    Z = copy(df)
    Z = select!(Z,  ["score_0", "score_1", "score_2"])
#    Z = Matrix(Z)

    # use Y = ["i_close_0", "i_close_1", "i_close_2"] to create a variable T
    Y = copy(df)
    Y = select!(Y, ["i_close_0", "i_close_1", "i_close_2"])

    # construct T variable
    N = length(Y[:, 1]) # number of observations
    Kx = length(X[1, :])  # number of time-invariant observations
    Kz = length(Z[1, :]) # number of time-varying observations

    T = zeros(length(Y[:, 1]))

    for n = 1:N
        if Y[n, 1] == 1
            T[n] = 1
        elseif Y[n, 1] == 0 && Y[n, 2] == 1
            T[n] = 2
        elseif Y[n, 1] == 0 && Y[n, 2] == 0 && Y[n, 3] == 1
            T[n] = 3
        elseif Y[n, 1] == 0 && Y[n, 2] == 0 && Y[n, 3] == 0
            T[n] = 4
        end
    end
    return X, Z, T, Y, N, Kx, Kz
end

X, Z, T, Y, N, Kx, Kz = prepare_data(; filename)

# res = Initialize(; N = N, K = K)

# Task 1. Quadrature method

# Nodes and weights for one dimensional integral
d1 = DataFrame(CSV.File(file_weigths_dim1))
nod1 = d1[:, "nodes"]
wt1 = d1[:, "weights"]

# Nodes and weights for two dimensional integral
d2 = DataFrame(CSV.File(file_weigths_dim2))
nod21 = d2[:, "nodes_1"]
nod22 = d2[:, "nodes_2"]
wt2 = d2[:, "weights"]

# transformation of nodes, so that they were for our region
function transform_nodes(; nodes1::Vector{Float64}, nodes21::Vector{Float64}, nodes22::Vector{Float64},  X::Vector{Real}, Z::Vector{Float32}, a0::Float64, a1::Float64, a2::Float64,
    β::Vector{Float64}, γ::Vector{Float64}, ρ::Float64)

    # question - should I create nodes for each individual separately???
    b = a0 .+ transpose(β)*X .+ transpose(γ)*Z
    nodes1_tr = log.(nodes1) .+ b
    nodes21_tr = log.(nodes21) .+ b

    b1 = a1 .+ transpose(β)*X .+ transpose(γ)*Z
    nodes22_tr = log.(nodes22) .+ b

    inv1 = 1 ./ nodes1
    inv21 = 1 ./ nodes21
    inv22 = 1 ./ nodes22

    return nodes1_tr, nodes21_tr, nodes22_tr, inv1, inv21, inv22
end



# Have to compute the nodes for each of individuals - for each individual.
x = [v for v in values(X[1000,:])]
z = [v for v in values(Z[1000,:])]

nodes1, nodes21, nodes22, inv1, inv21, inv22 = transform_nodes(; nodes1=nod1, nodes21=nod21, nodes22=nod22, X=x, Z=z, a0 = 0.0, a1 = -1.0, a2 = -1.0,
     β = zeros(Kx), γ = 0.3*ones(Kz), ρ= 0.5)


function compute_indiv_prob(; i::Int64, T::Float64, nodes1::Vector{Float64}, nodes21::Vector{Float64}, nodes22::Vector{Float64}, inv1::Vector{Float64}, inv21::Vector{Float64}, inv22::Vector{Float64},
    w1::Vector{Float64}, w2::Vector{Float64}, X::DataFrame, Z::DataFrame, a0::Float64, a1::Float64, a2::Float64,
    β::Vector{Float64}, γ::Vector{Float64}, ρ::Float64, σ0::Float64)

    # i - index stands for individual
    # T - duration 1-4

    dist = Normal(0, 1)

    x = [v for v in values(X[i,:])]
    z = [v for v in values(Z[i,:])]

    if T == 1.0
        P = cdf(dist, (-a0 - transpose(β)*x - transpose(γ)*z)/σ0)
    elseif T == 2.0
        n = length(nodes1)
        P = 0
        for j = 1:n
            ϵ0 = nodes1[j]
            P += cdf(dist, -a1 - transpose(β)*x - transpose(γ)*z - ρ*ϵ0)*pdf(dist, ϵ0/σ0)/σ0*inv1[j]*w1[j]
        end
    elseif T == 3.0 || T == 4.0
        n = length(nodes21)
        P = 0
        for j = 1:n
            ϵ0 = nodes21[j]
            ϵ1 = nodes22[j]
            if T == 3.0
                argum = -a2 - transpose(β)*x - transpose(γ)*z - ρ*ϵ1
            elseif T == 4.0
                argum = a2 + transpose(β)*x + transpose(γ)*z - ρ*ϵ1
            end
            P += cdf(dist, argum)*pdf(dist, ϵ1 - ρ*ϵ0)*pdf(dist, ϵ0/σ0)/σ0*inv21[j]*inv22[j]*w2[j]
        end
    end
    return P
end

ρ = 0.5
σ0 =  sqrt(1/(1-ρ^2))
### ADD RHO, SIGMA in parameters structure?


ind = 1000
dur = T[ind]
compute_indiv_prob(; i=ind, T = dur, nodes1 = nodes1, nodes21 = nodes21, nodes22 = nodes22, inv1 = inv1, inv21 = inv21, inv22 = inv22,
    w1 = wt1, w2 = wt2, X = X, Z = Z, a0 = 0.0, a1 = -1.0, a2 = -1.0, β = zeros(Kx), γ = 0.3*ones(Kz), ρ = 0.5, σ0 = 1.1547005383792515)

function compute_LL(; N::Int64, nodes1::Vector{Float64}, nodes21::Vector{Float64}, nodes22::Vector{Float64}, inv1::Vector{Float64}, inv21::Vector{Float64}, inv22::Vector{Float64},
    w1::Vector{Float64}, w2::Vector{Float64}, X::DataFrame, Z::DataFrame, T::Vector{Float64}, a0::Float64, a1::Float64, a2::Float64,
    β::Vector{Float64}, γ::Vector{Float64}, ρ::Float64, σ0::Float64)

    LL = 0
    P_ind = zeros(N)

    for i = 1:N
            ind = i
            dur = T[i]

            p_ind = compute_indiv_prob(; i=ind, T = dur, nodes1, nodes21, nodes22, inv1, inv21,inv22,
                w1, w2, X, Z, a0, a1, a2, β, γ, ρ, σ0)
            P_ind[i] = p_ind
                #    println(", p_ind = $p_ind\n")

            LL += log(p_ind)
    end

    return LL, P_ind
end

a0 = 0.0
a1 = -1.0
a2 = -1.0
β = zeros(Kx)
γ = 0.3*ones(Kz)
ρ = 0.5
σ0 = 1.1547005383792515

inv1
LL, P_ind = compute_LL(; N=N, nodes1=nodes1, nodes21=nodes21, nodes22=nodes22, inv1=inv1, inv21=inv21, inv22=inv22,
    w1=wt1, w2=wt2, X=X, Z=Z, T = T, a0 = a0, a1 = a1, a2 = a2,
        β = β, γ = γ, ρ = ρ, σ0 = σ0)

type1 = findall(x -> x == 1, T)
type2 = findall(x -> x == 2, T)
type3 = findall(x -> x == 3, T)
type4 = findall(x -> x == 4, T)

mean(P_ind[type1])
mean(P_ind[type2])
mean(P_ind[type3])
mean(P_ind[type4])

mean(P_ind[type1]) + mean(P_ind[type2]) + mean(P_ind[type3]) + mean(P_ind[type4])


dist = Normal(0, 1)
cdf(dist, 1.1) # use for probability Φ
pdf(dist, 0) # use for density ϕ


ρ = 0.5
σ0 = sqrt(1/(1-ρ)^2)
dist = Normal(0, 1)
pdf(dist, nodes1_tr[1]/σ0)/σ0
