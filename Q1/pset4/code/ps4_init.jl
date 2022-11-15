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

    σ::Float64   = 2.0  # coefficient of relative risk aversion
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

    T::Int64 = 50
    t_noss::Int64 = 30
end

mutable struct StationaryResults
    # value and policy functions, distribution are three-dimensional objects (assets - productivity - age)
    val_func::Array{Any, 3} # value function
    pol_func::Array{Any, 3} # policy function
    Γ::Array{Any, 3}        # distribution
    labor::Array{Any, 3}    # optimal labor choice

    # endogenous prices
    w::Float64 # wage
    r::Float64 # interest rate
    b::Float64 # pension benefits

    # aggregates
    K::Float64 # aggregate capital
    L::Float64 # aggregate labor
end

mutable struct TransitionResults
    # transition path of value and policy functions, distribution, and labor are four-dimensional objects (assets - productivity - age - time)
    val_func::Array{Float64,4}  # transition path of value function
    pol_func::Array{Float64,4}  # transition path of policy function
    Γ::Array{Float64,4}         # transition path of distribution
    labor::Array{Float64,4}     # transition path of labor    

    # transition path of prices
    w::Array{Float64,1} # transition path of wage
    r::Array{Float64,1} # transition path of interest rate
    b::Array{Float64,1} # transition path of SS pension benefits

    # transition path of aggregates
    K::Array{Float64} # transition path of aggregate capital
    L::Array{Float64} # transition path of aggregate labor

    # transition path of proportional labor income tax
    θ_path::Array{Float64}
end

# Initialize model primitives and stationary results
function Initialize(al::Float64, au::Float64, na::Int64, θ::Float64, z::Vector{Float64}, γ::Float64, T::Int64, t_noss::Int64)
    prim = primitives(al = al, au = au, na = na, θ = θ, z = z, γ = γ, T = T, t_noss = t_noss) #initialize primitives

    val_func = zeros(prim.na, 2, prim.N) # initialize value function
    pol_func = zeros(prim.na, 2, prim.N) # initialize policy function
    Γ        = zeros(prim.na, 2, prim.N) # initialize distribution
    labor    = zeros(prim.na, 2, prim.N) # initialize labor

    w = 1.05 
    r = 0.05
    b = 0.2
    if prim.θ == 0.0
        b = 0.0
    end

    K = 0.0
    L = 0.0

    res = StationaryResults(val_func, pol_func, Γ, labor, w, r, b, K, L) # initialize results struct

    prim, res
end

# Initialize model primitives and transition path results
function InitializeTrans(al::Float64, au::Float64, na::Int64, θ::Float64, z::Vector{Float64}, γ::Float64, T::Int64, t_noss::Int64)
    prim = primitives(al = al, au = au, na = na, θ = θ, z = z, γ = γ, T = T, t_noss = t_noss) #initialize primitives

    val_func = zeros(prim.T, prim.na, 2, prim.N) # initialize value function transition path
    pol_func = zeros(prim.T, prim.na, 2, prim.N) # initialize policy function transition path
    Γ        = zeros(prim.T, prim.na, 2, prim.N) # initialize distribution transition path
    labor    = zeros(prim.T, prim.na, 2, prim.N) # initialize labor transition path

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

    res = TransitionResults(val_func, pol_func, Γ, labor, w, r, b, K, L, θ_path) #initialize results struct

    prim, res
end
