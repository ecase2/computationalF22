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
    γ::Float64 = 0.42   # weight on consumption

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

function initInputs(θ_input::Float64, T_input::Int64, al::Int64, au::Int64, na_input::Int64, t_noss::Int64)
    θ      = θ_input
    T      = T_input
    t_noss = t_noss
    na = na_input
    a_grid = collect(range(au, al, na))

    ins    = inputs(θ, T, t_noss, na, a_grid)
    return ins
end

mutable struct results
    # results of the model
    # value and policy functions, distribution are four-dimensional objects (assets - productivity - age - period)
    val_func::Array{Float64, 4} # value function
    pol_func::Array{Float64, 4} # policy function
    Γ::Array{Float64, 4}        # distribution
    labor::Array{Float64, 4}    # optimal labor choice

    # endogenous prices, rows = t
    w::Array{Float64, 1} # wage transition path
    r::Array{Float64, 1} # interest rate
    b::Array{Float64, 1} # pension benefits

    # aggregates
    K::Array{Float64, 1} # aggregate capital
    L::Array{Float64, 1} # aggregate labor
end

include("ps4_ss_values.jl")
function initResults(prim::primitives, ins::inputs)
    @unpack na, T = ins
    @unpack N     = prim

    w = zeros(T)
    w[1]   = 1.335  # steady state wage with SS
    w[end] = 1.436  # steady state wage with no SS

    r = zeros(T)
    r[1]   = 0.037  # steady state interest rate with SS
    r[end] = 0.026  # steady state interest rate with no SS

    b = zeros(T)
    b[1] = 0.231    # steady state benefits with SS

    K = zeros(T)
    K[1]   = 2.953  # steady state aggregate capital with SS
    K[end] = 3.853  # steady state aggregate capital with no SS

    L = zeros(T)
    L[1]   = 0.383  # steady state aggregate labor with SS
    L[end] = 0.408  # steady state aggregate labor with no SS

    val_func = zeros(T, na, 2, N)
    pol_func = zeros(T, na, 2, N)
    labor    = zeros(T, na, 2, N)
    Γ = zeros(T, na, 2, N)

    val_worker, pol_worker, labor[1,:,:,:] = bellman_worker(prim, ins, w[1], r[1], 0.11)
    val_func[1, :, :, :], pol_func[1, :, :, :] = bellman_retiree(prim, ins, w[1], r[1], b[1], 0.11, val_worker, pol_worker)

    Γ[1, :,:,:]   = get_ss_distr(pol_func[1, :, :, :], prim, ins)

    val_worker, pol_worker, labor[end,:,:,:] = bellman_worker(prim, ins, w[end], r[end], 0.0)
    val_func[end, :, :, :], pol_func[end, :, :, :] = bellman_retiree(prim, ins, w[end], r[end], b[end], 0.0, val_worker, pol_worker)
    Γ[end, :,:,:] = get_ss_distr(pol_func[end, :, :, :],  prim, ins)

    res = results(val_func, pol_func, Γ, labor, w, r, b, K, L)
    return res
end

### For exercise 1 use t_noss = 2, while for exercise 2 use t_noss = 21
function initialize(θ_input::Float64 = 0.011, T_input::Int64 = 30, al::Int64 = 0, au::Int64 = 80, na_input::Int64 = 500, t_noss::Int64 = 1)
    prim = primitives()
    ins  = initInputs(θ_input, T_input, al, au, na_input, t_noss)
    res  = initResults(prim, ins)
    return prim, ins, res
end
