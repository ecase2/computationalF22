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
    # values we will not need to change
    N::Int64   = 66     # age of death
    n::Float64 = 0.011  # population growth
    R::Int64   = 46     # age of retirement

    σ::Int64   = 2      # coefficient of relative risk aversion
    α::Float64 = 0.36   # capital share
    δ::Float64 = 0.06   # depreciation rate
    β::Float64 = 0.97   # discount rate
    γ::Float64 = 0.42             # weight on consumption

    # transition matrix
    π_HH::Float64 = 0.9261
    π_LL::Float64 = 0.9811
    π::Matrix{Float64} = [π_HH (1-π_HH) ; (1-π_LL) π_LL] # productivity persistence probability matrix
    π0::Matrix{Float64} = [0.2037 0.7963]                # ergodic productivity distribution
    η::Array{Float64,1} = DataFrame(CSV.File(root*"/ef.csv"))[:, 1] # age-efficiency profile

    # normalizing the μ distribution - initial value
    μ_1::Float64 = (n*(1+n)^(N-1)) / ((1+n)^N - 1)

    z::Vector{Float64} = [3.0, 0.5]         # productivity shocks
    e::Matrix{Float64} = η*transpose(z)     # productivity
end

mutable struct inputs   # EMILY: is "changeables" or some other name more appropriate?
    # parameters or inputs that will change
    θ::Float64              # proportional labor income tax (for social security)
    T::Int64                # number of time periods
    t_noss::Int64           # time period when no SS takes effect

    # assets grid
    na::Int64   # number of points in the asset grid
    a_grid::Array{Float64} # asset grid
end

function initInputs(θ_input::Float64, T_input::Int64, al::Int64, au::Int64, na_input::Int64)
    a_grid = collect(range(au, al, na_input))
    t_noss = 21
    ins    = inputs(θ_input, T_input, t_noss, na_input, a_grid)
    return ins
end

mutable struct results
    # results of the model
    # value and policy functions, distribution are three-dimensional objects (assets - productivity - age)
    val_func::Array{Float64, 3} # value function
    pol_func::Array{Float64, 3} # policy function
    Γ::Array{Float64, 3}        # distribution
    labor::Array{Float64, 3}    # optimal labor choice

    # endogenous prices, rows = t
    w::Array{Float64, 1} # wage transition path
    r::Array{Float64, 1} # interest rate
    b::Array{Float64, 1} # pension benefits

    # aggregates
    K::Array{Float64, 1} # aggregate capital
    L::Array{Float64, 1} # aggregate labor
end

function initResults(prim::primitives, ins::inputs)
    @unpack na, T = ins
    @unpack N     = prim
    val_func = zeros(na, 2, N)
    pol_func = zeros(na, 2, N)
    Γ        = zeros(na, 2, N)
    labor    = zeros(na, 2, N)

    w = zeros(T)
    r = zeros(T)
    b = zeros(T)

    K = zeros(T)
    L = zeros(T)
    res = results(val_func, pol_func, Γ, labor, w, r, b, K, L)
    return res
end

function initialize(θ_input::Float64, T_input::Int64 ; al::Int64 = 0, au::Int64 = 80, na_input::Int64 = 500)
    prim = primitives()
    ins  = initInputs(θ_input, T_input, al, au, na_input)
    res  = initResults(prim, ins)
    return prim, ins, res
end
