
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
    𝐯::Array{Float64, 2} # value function
    𝐠::Array{Float64, 2} # decision rules
    𝛍::Array{Float64, 2} # cross-sectional distributions 
end

function initialize()
    θ = parameters()
    val_func      = zeros(θ.na, 2)
    decision_func = zeros(θ.na, 2)
    dist_func     = zeros(θ.na, 2)
    res = results(1.0, val_func, decision_func, dist_func)
    return θ, res
end

function utility(c::Float64; param::parameters = θ)
    @unpack 𝛂 = param
    (c^(1-𝛂) - 1)/(1 - 𝛂)
end

function bellman(res::results; param::parameters = θ)
    @unpack qstar, 𝐯, 𝐠 = res 
    @unpack 𝐒, 𝐀, na, a̲, 𝛃, 𝚷 = param
    vp = zeros(na, 2)   #value tomorrow (v_next)

    for si = 1:2    #loop through employment states 
        s = 𝐒[si]

        for ai = 1:na   #loop through today's a value 
            a = 𝐀[ai]
            candidate_max = -Inf #bad candidate max

            for api = 1:na  #loop through potential choices of tomorrow's a 
                ap = 𝐀[api]
                c = s+a - qstar*ap 
                if c >= 0    #check for positive consumption
                    val = utility(c) + 𝛃*𝚷[si, 1]*vp[api, 1] + 𝛃*𝚷[si, 2]*vp[api, 2]
                    if val > candidate_max
                        candidate_max = val #update the v amount to beat! 
                        𝐠[ai, si]     = ap  #update the corresponding choice of tmrw's a
                    end
                end
            end
            # now we know the best choice given today's s and a, so we should update the vp function 
            vp[ai, si] = candidate_max 
        end
    end
    return vp 
end

function v_iterate(res::results ;param::parameters = θ)
    @unpack ε = param 

    n   = 0     #count interations 
    err = 100   #initialize error 

    while err > ε
        vp = bellman(res)
        err = abs.(maximum(vp[:,1] .- res.𝐯[:,1] )) / abs.(maximum(vp[:,1])) + abs.(maximum(vp[:,2] .- res.𝐯[:,2] )) / abs.(maximum(vp[:,2]))
        res.𝐯 = vp 
        n += 1
        println("iteration $n with error $err")
    end
    println("THAT MF CONVERGED!!!! it took $n iterations tho")
end

function Tstar(res::results)
    @unpack 𝐠 = res 
end

function goodmarket(res::results; param::parameters = θ)
    @unpack 𝛍, 𝐠, qstar = res
    @unpack 𝐀, 𝐒 na = param

    sum = 0 
    for si = 1:2
        s = 𝐒[si] 
        for ai =1:na 
            # a  = 𝐀[ai]
            # ap = 𝐠[ai, si]
            # c  = s+a - qstar*ap 
            # sum += (c - s)*𝛍[ai, si]
            sum += 𝐠[ai,si]*𝛍[ai, si]
        end
    end
    return sum 
end

function assetmarket(res::results; param::parameters = θ)
    @unpack 𝛍, 𝐠 = res
    @unpack 𝐀, 𝐒 na = param

    sum = 0 
    for si = 1:2
        s = 𝐒[si]
        for ai = 1:na
            sum += 𝐠[ai, si]*𝛍[ai,si]
        end
    end
    return sum
end

function solve_model(res::results)
    ed = 100 # excess demand
    while abs(ed) > 1e-4
        v_iterate(res)
        #! something here about calculating the mu's 
        ed = assetmarket(res) #update excess demand amount
        if goodmarket(res) >0 #update guess for q 
            q = res.qstar + (1-res.qstar)/2 
            res.qstar = q
        elseif goodmarket(res) <0
            q = res.qstar - (1-res.qstar)/2 
            res.qstar = q
        end
    end
end