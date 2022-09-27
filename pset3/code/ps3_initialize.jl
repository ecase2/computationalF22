#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 3
    AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    CONTENTS:    This file initializes parameters and objects.
=#
@with_kw struct parameters
    N::Int64   = 66     # age of death
    n::Float64 = 0.011  # population growth
    R::Int64   = 46     # age of retirement

    σ::Int64   = 2      # coefficient of relative risk aversion
    α::Float64 = 0.36   # capital share
    δ::Float64 = 0.06   # depreciation rate
    β::Float64 = 0.97   # discount rate

    nz::Int64 = 2       # number of productivity shocks    
    η::Array{Float64,1} = DataFrame(CSV.File(root*"/ef.csv"))[:, 1] # age-efficiency profile

    # transition matrix
    π_HH::Float64 = 0.9261
    π_LL::Float64 = 0.9811
    π::Matrix{Float64} = [π_HH (1-π_HH) ; (1-π_LL) π_LL] # productivity persistence probability matrix
    π0::Matrix{Float64} = [0.2037 0.7963] 

    # normalizing the μ distribution - initial value
    μ_1::Float64 = (n*(1+n)^(N-1)) / ((1+n)^N - 1)

    # assets grid
    al::Int64 = 0 # lower bound for assets
    au::Int64 = 50 # upper bound for assets - some arbitraty number
    na::Int64 = 400 # number of points in the asset grid
    step::Float64 = (au - al)/(na-1)
    a_grid::Array{Float64} = collect(al:step:au)

    # convergence speed 
    λ::Float64 = 0.5
end

mutable struct results
    # I think that initializing using Any may slow down the process
    # parameters that change based on policy experiment
    θ::Float64              # proportional labor income tax
    z::Array{Float64, 1}    # productivity shocks
    e::Array{Float64, 1}    # productivity
    γ::Float64              # weight on consumption

    # value and policy functions, distribution are three-dimensional objects (assets - productivity - age)
    val_func::Array{Any, 3} # value function
    pol_func::Array{Any, 3} # policy function
    F::Array{Any, 3}        # distribution
    l::Array{Any, 3}             # optimal labor choice

    # the dimensions are [a, z, age, a']:
    c_grid::Array{Any, 4}        # consumption choices grid
    l_grid::Array{Any, 4}        # labor choices grid
       
    # endogenous prices
    w::Float64 # wage
    r::Float64 # interest rate
    b::Float64 # pension benefits

    # aggregates
    K::Float64 # aggregate capital
    L::Float64 # aggregate labor
end

#function for initializing model primitives and results
function Initialize(θ::Float64, z::Array{Float64, 1}, γ::Float64)
    par = parameters() #initialize primitives

    e = η*transpose(z)

    val_func = zeros(par.na, par.nz, par.N) #initial value function guess
    pol_func = zeros(par.na, par.nz, par.N) #initial policy function guess
    F = zeros(par.na, par.nz, par.N)
    c_grid = zeros(par.na, par.nz, par.N, par.na)
    l_grid = zeros(par.na, par.nz, par.N, par.na)
    l = zeros(par.na, par.nz, par.N)

    w = 1.05
    r = 0.05
    b = 0.2

    K = 0.0
    L = 0.0
    
    res  = results(θ, γ, z, e, val_func, pol_func, F, c, l, w, r, b, K, L) #initialize results struct

    prim, res
end