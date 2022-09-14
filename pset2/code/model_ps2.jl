
@with_kw struct parameters      # main parameters are bolded 
    𝛃::Float64  = 0.9932        # discount factor \bfbeta
    𝛂::Float64  = 1.5           # risk aversion coeff \bfalpha
    e::Float64  = 1.0           # employed earnings 
    u::Float64  = 0.5           # unemployed earnings 
    𝐒::Vector{Float64} = [e, u] # set of possible earnings \bfS

    ee::Float64 = 0.97      # π(s'= e | s = e)
    uu::Float64 = 0.5       # π(s'= u | s = u)
    𝚷::Matrix{Float64} = [[ee, 1-uu] [1-ee, uu]]  # transition matrix (\bfPi)

    a̲::Float64  = -2        # asset holdings lower bound (a\underbar)
    ā::Float64  = 5         # asset holdings upper bound (a\bar)
    na::Int64 = 1000        # number of asset grid points
    𝐀::Array{Float64,1} = collect(range(a̲, length = na, stop = ā)) # A grid 

    ε::Float64  = 1e-4      # tolerance for convergence (\varepislon)
end


mutable struct results 
    qstar::Float64       # market clearing price ∈ [0,1]
    𝐠::Array{Float64, 2} # decision rules
    𝛍::Array{Float64, 2} # cross-sectional distributions 
end

function initialize()
    θ = parameters()


    results(0.5, zeros(θ.na,2), zeros(θ.na,2) )
end

function utility(c::Float64; param::parameters = θ)
    @unpack 𝛂 = param
    (c^(1-𝛂) - 1)/(1 - 𝛂)
end

function bellman()
    
end