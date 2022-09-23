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

    π_HH::Float64 = 0.9261
    π_LL::Float64 = 0.9811
    π::Matrix{Float64} = [π_HH (1-π_HH) ; (1-π_LL) π_LL] # productivity persistence probability matrix

    w::Float64 = 1.05   # wage
    r::Float64 = 0.05   # interest rate
    b::Float64 = 0.2    # pension benefits

    # normalizing the μ distribution - initial value 
    μ_1::Float64 = (n*(1+n)^(N-1)) / ((1+n)^N - 1)
end

mutable struct results
    # unsure of structure yet
end