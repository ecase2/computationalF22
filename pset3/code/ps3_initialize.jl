#=====================================
    INITIALIZE PARAMETERS AND OBJECTS
=====================================#
mutable struct parameters
    N::Int64    # age of death
    n::Float64  # population growth
    R::Int64    # age of retirement
    θ::Float64  # proportional labor income tax
    γ::Float64  # weight on consumption
    σ::Int64    # coefficient of relative risk aversion
    α::Float64  # capital share
    δ::Float64  # depreciation rate
    β::Float64  # discount rate

    zᴴ::Float64     # high productivity value
    zᴸ::Float64     # low productivity value
    z::Array{Float64}
    nz::Int64       # number of productivity shocks

    π_HH::Float64
    π_LL::Float64
    π::Matrix{Float64} # productivity persistence probability matrix

    # normalizing the μ distribution - initial value
    μ_1::Float64

    # assets grid
    al::Int64   # lower bound for assets
    au::Int64   # upper bound for assets - some arbitraty number
    na::Int64   # number of points in the asset grid
    step::Float64
    a_grid::Array{Float64}
end
function param_init(modeltype::String)
    # the modeltype == "benchmark" values are:
    N   = 66     # age of death
    n   = 0.011  # population growth
    R   = 46     # age of retirement
    θ   = 0.11   # proportional labor income tax
    γ   = 0.42   # weight on consumption
    if modeltype == "exogenouslabor"
        γ = 1.0 
    end 
    σ   = 2      # coefficient of relative risk aversion
    α   = 0.36   # capital share
    δ   = 0.06   # depreciation rate
    β   = 0.97   # discount rate

    zᴴ  = 3.0    # high productivity value
    zᴸ  = 0.5    # low productivity value
    if modeltype == "norisk"
        zᴴ = 0.5 
        zᴸ = 0.5
    end
    z   = [zᴴ, zᴸ]
    nz  = 2      # number of productivity shocks

    π_HH = 0.9261
    π_LL = 0.9811
    π = [π_HH (1-π_HH) ; (1-π_LL) π_LL] # productivity persistence probability matrix

    # normalizing the μ distribution - initial value
    μ_1 = (n*(1+n)^(N-1)) / ((1+n)^N - 1)

    # assets grid
    al = 0 # lower bound for assets
    au = 50 # upper bound for assets - some arbitraty number
    na = 400 # number of points in the asset grid
    step = (au - al)/(na-1)
    a_grid = collect(al:step:au)

    par = parameters(N,n,R,θ,γ,σ,α,δ,β,zᴴ,zᴸ,z,nz,π_HH,π_LL,π,μ_1,al,au,na,step,a_grid)
    return par
end

mutable struct results
    # unsure of structure yet
    # I think that initializing using Any may slow down the process
    # value and policy functions, distribution are three-dimensional objects (assets - productivity - age)
    val_func::Array{Any, 3} # value function
    pol_func::Array{Any, 3} # policy function
    F::Array{Any, 3}        # distribution

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
    # the dimensions are [a, z, age, a']:
    c::Array{Any, 4}        # consumption choices
    l::Array{Any, 4}        # labor choices 
    # μ::Array{Float64, 1}    # cohort size by age starting with μ_1  
end

#function for initializing model primitives and results
function Initialize(modeltype::String = "benchmark")

    par = param_init(modeltype) #initialize primtiives

    val_func = zeros(par.na, par.nz, par.N) #initial value function guess
    pol_func = zeros(par.na, par.nz, par.N) #initial policy function guess
    F = zeros(par.na, par.nz, par.N)

    w = 1.05
    r = 0.05
    b = 0.2

    e = DataFrame(CSV.File(root*"/ef.csv"))[:, 1]*transpose(par.z)

    c = zeros(par.na, par.nz, par.N, par.na)
    l = zeros(par.na, par.nz, par.N, par.na)
    

    res = results(val_func, pol_func, F, w, r, b, e) #initialize results struct
    grid = grids(c, l)
    par, res, grid
end
