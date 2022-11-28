#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 2 (Q2)
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#


function prepare_data(; filename::String)
    df = DataFrame(load(filename))

    # Observed characteristics - time-invariant
    X = copy(df)
    X = select!(X,  ["score_0", "rate_spread", "i_large_loan", "i_medium_loan", "i_refinance", "age_r",  "cltv",
                "dti", "cu", "first_mort_r", "i_FHA",
                 "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"])

    # Add the column for constants
    X_add = ones(length(X[:, 1]))
    X = hcat(X, X_add)
    X = X .+ 0.0                            # to convert to type Float64

    # Observed characteristics - time-invariant
    Z = copy(df)
    Z = select!(Z,  ["score_0", "score_1", "score_2"])
    Z = Z .+ 0.0                            # to convert to type Float64

    Y = copy(df)
    Y = select!(Y, ["i_close_0", "i_close_1", "i_close_2"])

    # construct T variable
    N = length(Y[:, 1])                     # number of observations
    Kx = length(X[1, :])                    # number of time-invariant observations
    Kz = length(Z[1, :])                    # number of time-varying observations

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


@with_kw struct Params
    a0::Float64 = 0.0
    a1::Float64 = -1.0
    a2::Float64 = -1.0
    β::Vector{Float64} = zeros(Kx)
    γ::Vector{Float64} = 0.3*ones(Kz)
    ρ::Float64 = 0.5
    σ0::Float64 = sqrt(1/(1-ρ)^2)
end

function Initialize(; Kx::Int64, Kz::Int64)
    par = Params(β = zeros(Kx) , γ = 0.3*ones(Kz))

end


####################################### Quadrature algorithm ############################33


function read_nodes(;file_weigths_dim1::String, file_weigths_dim2::String)

    # Nodes and weights for one dimensional integral
    d1 = DataFrame(CSV.File(file_weigths_dim1))
    nod1 = d1[:, "nodes"]
    w1 = d1[:, "weights"]

    # Nodes and weights for two dimensional integral
    d2 = DataFrame(CSV.File(file_weigths_dim2))
    nod21 = d2[:, "nodes_1"]
    nod22 = d2[:, "nodes_2"]
    w2 = d2[:, "weights"]

    return nod1, nod21, nod22, w1, w2

end

# Transformation of nodes, so that they were for our region (for one individual)
function transform_nodes(par::Params; nodes1::Vector{Float64}, nodes21::Vector{Float64}, nodes22::Vector{Float64},  X::Vector{Float64}, Z::Vector{Float64},
    a0::Float64, a1::Float64, a2::Float64, β::Vector{Float64}, γ::Vector{Float64}, ρ::Float64)

    σ0 = sqrt(1/(1-ρ)^2)

    b = a0 .+ transpose(β)*X .+ γ[1]*Z[1]
    nod1_tr = log.(nodes1) .+ b
    nod21_tr = log.(nodes21) .+ b

    b1 = a1 .+ transpose(β)*X .+ γ[2]*Z[2]
    nod22_tr = log.(nodes22) .+ b1

    inv1 = 1 ./ nodes1
    inv21 = 1 ./ nodes21
    inv22 = 1 ./ nodes22

    return nod1_tr, nod21_tr, nod22_tr, inv1, inv21, inv22
end


function compute_indiv_prob_quadr(par::Params; i::Int64, T::Float64, nodes1::Vector{Float64}, nodes21::Vector{Float64}, nodes22::Vector{Float64}, inv1::Vector{Float64}, inv21::Vector{Float64}, inv22::Vector{Float64},
    w1::Vector{Float64}, w2::Vector{Float64}, x::Vector{Float64}, z::Vector{Float64}, a0::Float64,
    a1::Float64, a2::Float64, β::Vector{Float64}, γ::Vector{Float64}, ρ::Float64)

    σ0 = sqrt(1/(abs(1-ρ)^2))

    # i - index stands for individual
    # T - duration 1-4

    dist = Normal(0, 1)

    if T == 1.0
        P = cdf(dist, (-a0 - transpose(x)*β - z[1]*γ[1])/σ0)
    elseif T == 2.0
        n = length(nodes1)
        P = 0
        for j = 1:n
            ϵ0 = nodes1[j]
            P += cdf(dist, -a1 - transpose(x)*β - z[2]*γ[2] - ρ*ϵ0)*pdf(dist, ϵ0/σ0)/σ0*inv1[j]*w1[j]
        end
    elseif T == 3.0 || T == 4.0
        n = length(nodes21)
        P = 0
        for j = 1:n
            ϵ0 = nodes21[j]
            ϵ1 = nodes22[j]
            if T == 3.0
                argum = -a2 - transpose(x)*β -  z[3]*γ[3] - ρ*ϵ1
            elseif T == 4.0
                argum = a2 + transpose(x)*β +  z[3]*γ[3] - ρ*ϵ1
            end
            P += cdf(dist, argum)*pdf(dist, ϵ1 - ρ*ϵ0)*pdf(dist, ϵ0/σ0)/σ0*inv21[j]*inv22[j]*w2[j]
        end
    end
    return P
end


function quadr_LL(par::Params; N::Int64, nodes1::Vector{Float64}, nodes21::Vector{Float64}, nodes22::Vector{Float64},
    w1::Vector{Float64}, w2::Vector{Float64}, X::DataFrame, Z::DataFrame, T::Vector{Float64}, a0::Float64,
    a1::Float64, a2::Float64, β::Vector{Float64}, γ::Vector{Float64}, ρ::Float64)

    LL = 0
    P_ind = zeros(N)

    for i = 1:N
            ind = i
            dur = T[i]

            # transform the nodes for an individual
            x = [v for v in values(X[i,:])]
            z = [v for v in values(Z[i,:])]
            nod1_tr, nod21_tr, nod22_tr, inv1, inv21, inv22 = transform_nodes(par; nodes1, nodes21, nodes22, X = x, Z = z,
             a0 = a0, a1 = a1, a2 = a2, β = β, γ = γ, ρ = ρ)

            # find probability for a given individual
            p_ind = compute_indiv_prob_quadr(par; i=ind, T = dur, nodes1 = nod1_tr, nodes21 = nod21_tr, nodes22 = nod22_tr, inv1 = inv1, inv21 = inv21,
                inv22 = inv22, w1 = w1, w2 = w2, x =  x, z = z, a0 = a0, a1 = a1, a2 = a2, β = β, γ = γ, ρ = ρ)
            P_ind[i] = p_ind

            LL += log(p_ind)
    end

    return LL
end


################################### GHK algorithm #############################################333

### NOT SURE :((( NEED TO TRY DOING THE SAME BUT WITH HALTON SEQUENCE
function compute_indiv_prob_ghk(par::Params; i::Int64, T::Float64, x::Vector{Float64}, z::Vector{Float64}, R::Int64)
    @unpack a0, a1, a2, β, γ, ρ, σ0 = par

    distN = Normal(0, 1)
    distU = Uniform(0, 1)

    if T == 1.0
        P = cdf(distN, (-a0 - transpose(x)*β -  z[1]*γ[1])/σ0)
    elseif T == 2.0

        draws = rand(distU, R)
        CDF = cdf(distN, a0 + transpose(x)*β +  z[1]*γ[1])
        ϵ0_grid = σ0*quantile(distN, draws.*CDF)

        P = 0
        for j = 1:R
            ϵ0 = ϵ0_grid[j]
            P += cdf(distN, -a1 - transpose(x)*β -  z[2]*γ[2] - ρ*ϵ0)*CDF
        end
        P = 1/R*P

    elseif T == 3.0 || T == 4.0

        draws0 = rand(distU, R)
        draws1 = rand(distU, R)

        CDF0 = cdf(distN, a0 + transpose(x)*β +  z[1]*γ[1])
        ϵ0_grid = σ0*quantile(distN, draws0.*CDF0)

        CDF1 = cdf(distN, a1 .+ transpose(x)*β .+  z[2]*γ[2] .+ ρ*ϵ0_grid)
        ϵ1_grid = quantile(distN, draws1.*CDF1)

        P = 0
        for j = 1:R
            ϵ0 = ϵ0_grid[j]
            ϵ1 = ϵ1_grid[j]
            if T == 3.0
                argum = -a2 - transpose(x)*β -  z[3]*γ[3] - ρ^2*ϵ0 - ρ*ϵ1
            elseif T == 4.0
                argum = a2 + transpose(x)*β + z[3]*γ[3] + ρ^2*ϵ0 + ρ*ϵ1
            end
            P += cdf(distN, argum)*CDF0*CDF1[j]
        end
        P = 1/(R^2)*P
    end
    return P

end


function ghk_LL(par::Params; X::DataFrame, Z::DataFrame, T::Vector{Float64}, N::Int64, R::Int64)
    # R is number of draws of random variables for each
    @unpack a0, a1, a2, β, γ, ρ, σ0 = par

    LL = 0
    P_ind = zeros(N)

    for i = 1:N
    #    println("i = $i")
            ind = i

            # transform the nodes for an individual
            x = [v for v in values(X[i,:])]
            z = [v for v in values(Z[i,:])]

            # find probability for a given individual
            p_ind = compute_indiv_prob_ghk(par;  i = i, T = T[i], x = x, z = z, R = R)
            P_ind[i] = p_ind

            LL += log(p_ind)
    end

    return LL
end

############################## Accept/Reject ############################################
function Accept_reject(par::Params; X::DataFrame, Z::DataFrame, T::Vector{Float64},  N::Int64, R::Int64)

    @unpack a0, a1, a2, β, γ, ρ, σ0 = par

    distU = Uniform()
    distN = Normal(0, 1)

    Random.seed!(1234)

    draw0 = rand(distU, R)
    draw1 = hcat(rand(distU, R), rand(distU, R))
    draw2 = hcat(rand(distU, R), rand(distU, R), rand(distU, R))

    e0, η1, η2 = zeros(R), zeros(R), zeros(R)

    AR = 0
    Acc = 0
    for i = 1:N
        x = [v for v in values(X[i,:])]
        z = [v for v in values(Z[i,:])]
        A = zeros(R)
        if T[i] == 1.0
            cutoff = -a0 - transpose(x)*β - z[1]*γ[1]
            for j = 1:R
                e0[j] = σ0*quantile.(distN, draw0[j])
                if e0[j] < cutoff
                    A[j] = 1
                end
            end
            Acc = mean(A)
        elseif T[i] == 2.0
            e0 = σ0*quantile.(distN, draw1[:, 1])
            η1 = quantile.(distN, draw1[:, 2])

            cutoff0 = a0 + transpose(x)*β + z[1]*γ[1]
            A0 = zeros(R)
            A1 = zeros(R)
            for j = 1:R
                if e0[j] < cutoff0
                    A0[j] = 1
                end
                cutoff1 = -a1 - transpose(x)*β - z[2]*γ[2] - ρ*e0[j]
                if η1[j] < cutoff1
                    A1[j] = 1
                end
            end
            Acc = mean(A0.*A1)
        elseif T[i] == 3.0 || T == 4.0
            e0 = σ0*quantile.(distN, draw2[:, 1])
            η1 = quantile.(distN, draw2[:, 2])
            η2 = quantile.(distN, draw2[:, 3])

            cutoff0 = a0 + transpose(x)*β + z[1]*γ[1]
            A0 = zeros(R)
            A1 = zeros(R)
            A2 = zeros(R)
            for j = 1:R
                cutoff1 = a1 + transpose(x)*β + z[2]*γ[2] - ρ*e0[j]

                if T[i] == 3.0
                    cutoff2 = -a2 - transpose(x)*β - z[3]*γ[3] - ρ^2*e0[j] - ρ*η1[j]
                elseif T[i] == 4.0
                    cutoff2 = a2 + transpose(x)*β + z[3]*γ[3] - ρ^2*e0[j] - ρ*η1[j]
                end

                if e0[j] < cutoff0
                    A0[j] = 1
                end
                if η1[j] < cutoff1
                    A1[j] = 1
                end
                if η2[j] < cutoff2
                    A2[j] = 1
                end
            end
            Acc = mean(A0.*A1.*A2)
        end
        AR = AR + log(Acc)
        if AR == -Inf
            AR = -10e-15
        end
    end
    return AR
end


#################################### For estimation
function loglike(params)

   a0 = params[1]
   a1 = params[2]
   a2 = params[3]
   β = params[4:19]
   γ = params[20:22]
   ρ = params[23]

   dist = Normal(0, 1)

   LL = 0

   for i = 1:N
           ind = i
           dur = T[i]

           x = [v for v in values(X[i,:])]
           z = [v for v in values(Z[i,:])]
           nod1_tr, nod21_tr, nod22_tr, inv1, inv21, inv22 = transform_nodes(par; nodes1 = nod1, nodes21 = nod21, nodes22 = nod22, X = x, Z = z,
           a0 = a0, a1 = a1, a2 = a2, β = β, γ = γ, ρ = ρ)

           p_ind = compute_indiv_prob_quadr(par; i=ind, T = dur, nodes1 = nod1_tr, nodes21 = nod21_tr, nodes22 = nod22_tr, inv1 = inv1, inv21 = inv21,
               inv22 = inv22, w1 = w1, w2 = w2, x =  x, z = z, a0 = a0, a1 = a1, a2 = a2, β = β, γ = γ, ρ = ρ)

           LL += log(p_ind)
   end
   return -LL
end
