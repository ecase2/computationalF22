#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 6
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

using Parameters, DataFrames

### Structure of model primitives (unmutable objects)
@with_kw struct Params
    β::Float64          = 0.8
    θ::Float64          = 0.64
    s::Array{Float64}   = [3.98e-4, 3.58, 6.82, 12.18, 18.79]
    Ns::Int64           = length(s)
    F::Array{Float64,2} = [0.6598 0.2600 0.0416 0.0331 0.0055;
                           0.1997 0.7201 0.0420 0.0326 0.0056;
                           0.2000 0.2000 0.5555 0.0344 0.0101;
                           0.2000 0.2000 0.2502 0.3397 0.0101;
                           0.2000 0.2000 0.2500 0.3400 0.0100]

    ν::Array{Float64}   = [0.37, 0.4631, 0.1102, 0.0504, 0.0063]
    A::Float64          = 1/200
    cf::Int64           = 10
    ce::Int64           = 5
end

mutable struct Results
    p::Float64  # price level
    I::Float64  # mass of incumbents
    M::Float64  # mass of entrants
    D::Float64  # mass of exits
    L::Float64  # aggregate labor
    L_I::Float64 # labor of incumbents
    L_E::Float64 # labor of entrants
    fLE::Float64 # fraction of labor in entrants
    μ::Array{Float64} # distribution of firms
    pol_func::Array{Float64} # exit policy function
    W::Array{Float64}
end

# Initialize
function Initialize(;cf::Int64)
    par = Params(cf = cf)

    p = 0.5 # price level
    I = 0.7  # mass of incumbents
    M = 0.3 # mass of entrants
    D = 0.1 # mass of exits
    L = 1   # aggregate labor
    L_I = 0.6 # labor of incumbents
    L_E = 0.4 # labor of entrants
    fLE = 0.2 # fraction of labor in entrants - WHAT IS IT

    μ = zeros(par.Ns) # distribution of firms
    pol_func =  zeros(par.Ns)
    W = zeros(par.Ns)

    res = Results(p, I, M, D, L, L_I, L_E, fLE, μ, pol_func, W)   # initialize results structure

    return par, res
end


### Firm's problem

function firms_problem(par::Params, res::Results; price::Float64)
    @unpack s, Ns, F, β, θ, cf = par
    @unpack W, pol_func = res

    W_temp = zeros(Ns)
    exit = zeros(Ns)

    for s_index = 1:Ns
        n = (θ*price*s[s_index])^(1/(1-θ))
        prof = θ*s[s_index]*n^θ - n - price*cf


        W_future = β*sum(transpose(F[s_index, :])*W_temp)

        if W_future < 0
            W_temp[s_index] = prof
            exit[s_index] = 1
        else
            W_temp[s_index] = prof + W_future
            exit[s_index] = 0
        end
    end
    res.pol_func = exit

    return W_temp
end


#function to iterate over value function
function v_iterate(par::Params, res::Results; tol::Float64=1e-3, price::Float64 = 0.7)
    #@unpack W = res
    err = 100.0

    while err > tol
        W_temp = firms_problem(par, res; price)
        err = maximum(abs.(res.W .- W_temp))
        res.W = W_temp
    end
end


### To check how it works
# par, res = Initialize(; cf = 10)
# @time v_iterate(par, res; price = 0.7)
