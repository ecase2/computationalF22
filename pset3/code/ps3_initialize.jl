#=====================================
    INITIALIZE PARAMETERS AND OBJECTS
=====================================#
@with_kw struct parameters
    N::Int64   = 66     # age of death
    n::Float64 = 0.011  # population growth
    R::Int64   = 46     # age of retirement
    θ::Float64 = 0.11   # proportional labor income tax
    γ::Float64 = 0.42   # weight on consumption
    σ::Int64   = 2      # coefficient of relative risk aversion
    α::Float64 = 0.36   # capital share
    δ::Float64 = 0.06   # depreciation rate
    β::Float64 = 0.97   # discount rate

    zᴴ::Float64 = 3.0   # high productivity value
    zᴸ::Float64 = 0.5   # low productivity value
    z::Array{Float64} = [zᴴ, zᴸ]
    nz::Int64 = 2 # number of productivity shocks

    pzᴴ::Float64 = 0.2037 # probability to born with high productivity
    pzᴸ::Float64 = 0.7963 # probability to born with low productivity

    π_HH::Float64 = 0.9261
    π_LL::Float64 = 0.9811
    π::Matrix{Float64} = [π_HH (1-π_HH) ; (1-π_LL) π_LL] # productivity persistence probability matrix

    #w::Float64 = 1.05   # wage
    #r::Float64 = 0.05   # interest rate
    #b::Float64 = 0.2    # pension benefits

    # normalizing the μ distribution - initial value
    μ_1::Float64 = (n*(1+n)^(N-1)) / ((1+n)^N - 1)

    # assets grid
    al::Int64 = 0 # lower bound for assets
    au::Int64 = 50 # upper bound for assets - some arbitraty number
    na::Int64 = 400 # number of points in the asset grid
    step::Float64 = (au - al)/(na-1)
    a_grid::Array{Float64} = collect(al:step:au)
end

mutable struct results
    # unsure of structure yet
    # I think that initializing using Any may slow down the process
    # value and policy functions, distribution are three-dimensional objects (assets - productivity - age)
    val_func::Array{Any, 3} # value function
    pol_func::Array{Any, 3} # policy function
    μ::Array{Any, 3} # distribution

    # endogenous prices
    w::Float64 # wage
    r::Float64 # interest rate

    # social security benefit
    b::Float64 # pension benefits

    # array to hold worker's productivity
    e::Array{Float64, 2} # dependendent on the deterministic age-efficiency profile and shock z

end

# structure to hold the grid for consumption and labor choices
# it should make the computation faster
# consumption and labor depends on current level of assets,  productivity, age, and tomorrow's level of assets
mutable struct grids
        c::Array{Any, 4}
        l::Array{Any, 4}
    #    u::Array{Any, 4}
end

#function for initializing model primitives and results
function Initialize()
    par = parameters() #initialize primtiives
    # NOT SURE ABOUT INITIALIZATION of value and policy function
    val_func = zeros(par.na, par.nz, par.N) #initial value function guess
    pol_func = zeros(par.na, par.nz, par.N) #initial policy function guess
    μ = zeros(par.na, par.nz, par.N)

    w = 1.05
    r = 0.05
    b = 0.2

    e = DataFrame(CSV.File(root*"/ef.csv"))[:, 1]*transpose(par.z)

    c = zeros(par.na, par.nz, par.N, par.na)
    l = zeros(par.na, par.nz, par.N, par.na)
    #u = zeros(par.na, par.nz, par.N, par.na)

    res = results(val_func, pol_func, μ, w, r, b, e) #initialize results struct
    grid = grids(c, l)
    par, res, grid
end
