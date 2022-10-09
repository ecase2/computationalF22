#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 4
    * AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    * CONTENTS:    Creates functions to initialize
    * NOTE:
    *    primitives - immutable struct, everything we don't need to change
    *    parameters - mutable struct, inputs we need to change
    *    results    - mutable struct, contains results of our model
=#
@with_kw struct primitives
    N::Int64   = 66     # age of death
    n::Float64 = 0.011  # population growth
    R::Int64   = 46     # age of retirement

    σ::Float64   = 2.0      # coefficient of relative risk aversion
    α::Float64 = 0.36   # capital share
    δ::Float64 = 0.06   # depreciation rate
    β::Float64 = 0.97   # discount rate

    η::Array{Float64,1} = DataFrame(CSV.File(root*"/ef.csv"))[:, 1] # age-efficiency profile

    # transition matrix
    π_HH::Float64 = 0.9261
    π_LL::Float64 = 0.9811
    π::Matrix{Float64} = [π_HH (1-π_HH) ; (1-π_LL) π_LL] # productivity persistence probability matrix
    π0::Matrix{Float64} = [0.2037 0.7963]

    # normalizing the μ distribution - initial value
    μ_1::Float64 = (n*(1+n)^(N-1)) / ((1+n)^N - 1)

    # assets grid
    al::Int64 = 0   # lower bound for assets
    au::Int64 = 20  # upper bound for assets - some arbitraty number
    na::Int64 = 200 # number of points in the asset grid
    step::Float64 = (au - al)/(na-1)
    a_grid::Array{Float64} = collect(al:step:au)

    θ::Float64 = 0.11
    z::Vector{Float64} = [3.0, 0.5]
    γ::Float64 = 0.42
    e::Array{Float64, 2} = η*transpose(z)

    T::Int64 = 1
    t_noss::Int64 = 2

end

mutable struct results

    # value and policy functions, distribution are three-dimensional objects (assets - productivity - age)
    val_func::Array{Any, 3} # value function
    pol_func::Array{Any, 3} # policy function
    Γ::Array{Any, 3}        # distribution
    labor::Array{Any, 3}        # optimal labor choice

    # endogenous prices
    w::Float64 # wage
    r::Float64 # interest rate
    b::Float64 # pension benefits

    # aggregates
    K::Float64 # aggregate capital
    L::Float64 # aggregate labor
end

#function for initializing model primitives and results
function Initialize(al::Float64, au::Float64, na::Int64, θ::Float64, z::Vector{Float64}, γ::Float64, T::Int64, t_noss::Int64)
    prim = primitives(al = al, au = au, na = na, θ = θ, z = z, γ = γ, T = T, t_noss = t_noss) #initialize primitives

    val_func = zeros(prim.na, 2, prim.N) #initial value function guess
    pol_func = zeros(prim.na, 2, prim.N) #initial policy function guess
    Γ        = zeros(prim.na, 2, prim.N)
    labor    = zeros(prim.na, 2, prim.N)

    w = 1.05
    r = 0.05
    b = 0.2
    if prim.θ == 0.0
        b = 0.0
    end

    K = 0.0
    L = 0.0

    res  = results(val_func, pol_func, Γ, labor, w, r, b, K, L) #initialize results struct

    prim, res
end


mutable struct results_trans
    # results of the model
    val_func::Array{Float64,4}  #value function is 4D - age, z1, z2, T
    pol_func::Array{Float64,4}
    labor::Array{Float64,4}  #labor supply function
    Γ::Array{Float64,4}

    # PRICES AND AGGS NEED TO BE ARRAYS NOW - THEY ARE PATHS OVER TIME
    w::Array{Float64,1} #wage path
    r::Array{Float64,1} #interest rate path
    b::Array{Float64,1} #SS pension benefits path
    K::Array{Float64} #agg capital path
    L::Array{Float64} #agg labor path

    #NEW
    θ_path::Array{Float64}
end


function Initialize_trans(al::Float64, au::Float64, na::Int64, θ::Float64, z::Vector{Float64}, γ::Float64, T::Int64, t_noss::Int64)
    prim = primitives(al = al, au = au, na = na, θ = θ, z = z, γ = γ, T = T, t_noss = t_noss) #initialize primitives

    val_func = zeros(prim.T, prim.na, 2, prim.N) #initial value function guess
    pol_func = zeros(prim.T, prim.na, 2, prim.N) #initial policy function guess
    Γ        = zeros(prim.T, prim.na, 2, prim.N)
    labor    = zeros(prim.T, prim.na, 2, prim.N)

    w = 1.05*ones(prim.T)
    r = 0.05*ones(prim.T)
    b = 0.2*ones(prim.T)

    K = zeros(prim.T)
    L = zeros(prim.T)

    θ_path = θ*ones(prim.T)

    for t = 1:T
        if t >= t_noss
            θ_path[t] = 0
            b[t] = 0
        end
    end

    res  = results_trans(val_func, pol_func, Γ, labor, w, r, b, K, L, θ_path) #initialize results struct

    prim, res
end
