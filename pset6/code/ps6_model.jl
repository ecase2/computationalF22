#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 6
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#
using Parameters
# Structure of model primitives (unmutable objects)
@with_kw struct Params
    β::Float64          = 0.8                                   # discount factor
    θ::Float64          = 0.64                                  # productivity parameter
    s::Array{Float64}   = [3.98e-4, 3.58, 6.82, 12.18, 18.79]   # productivity shocks
    Ns::Int64           = length(s)                             # number of shocks types
    F::Array{Float64,2} = [0.6598 0.2600 0.0416 0.0331 0.0055;  # transition matrix for shocks
                           0.1997 0.7201 0.0420 0.0326 0.0056;
                           0.2000 0.2000 0.5555 0.0344 0.0101;
                           0.2000 0.2000 0.2502 0.3397 0.0101;
                           0.2000 0.2000 0.2500 0.3400 0.0100]

    ν::Array{Float64}   = [0.37, 0.4631, 0.1102, 0.0504, 0.0063] # distribution of types for entrants
    A::Float64          = 1/200                                  # labor utility parameter
    cf::Int64           = 10                                     # costs of staying
    ce::Int64           = 5                                      # entry costs
    α::Int64            = 1                                      # EV1 distribution parameter
end

mutable struct Results
    p::Float64                                                   # price level
    m_inc::Float64                                               # mass of incumbents
    m_entr::Float64                                              # mass of entrants
    m_exit::Float64                                              # mass of exits
    L::Float64                                                   # aggregate labor
    L_inc::Float64                                               # labor of incumbents
    L_entr::Float64                                              # labor of entrants
    frac_L_entr::Float64                                         # fraction of labor in entrants
    μ::Array{Float64}                                            # distribution of firms
    pol_func::Array{Float64}                                     # exit policy function
    W::Array{Float64}                                            # value function for baseline case
    labor::Array{Float64}                                        # optimal labor demand

    ### Part for shocks
    U0::Array{Float64}                                           # ex-ante value (version with shocks)
    V0:: Array{Float64}                                          # value of staying
    V1::Array{Float64}                                           # value of exiting
    σ::Array{Float64}                                            # probability to exit = policy function in the version with shocks
end

# Function to initialize parameters and results structures
function Initialize(; cf::Int64, α::Int64)
    par = Params(cf = cf, α = α)

    p = 0.5
    m_inc = 0.7
    m_entr = 0.3
    m_exit = 0.1
    L = 1
    L_inc = 0.6
    L_entr = 0.4
    frac_L_entr = 0.2

    μ = zeros(par.Ns)
    pol_func =  zeros(par.Ns)
    W = zeros(par.Ns)
    labor = zeros(par.Ns)

    U0 = zeros(par.Ns)
    V0 = zeros(par.Ns)
    V1 = zeros(par.Ns)
    σ  = zeros(par.Ns)

    res = Results(p, m_inc, m_entr, m_exit, L, L_inc, L_entr, frac_L_entr, μ, pol_func, W, labor, U0, V0, V1, σ)

    return par, res
end


# Firm's problem for a baseline version

function firms_problem(par::Params, res::Results; price::Float64)
    @unpack s, Ns, F, β, θ, cf = par
    @unpack W, pol_func, labor = res

    W_temp = zeros(Ns)
    exit = zeros(Ns)

    for s_index = 1:Ns
        n = (1/(θ*price*s[s_index]))^(1/(θ-1))
        labor[s_index] = n

        prof = price*s[s_index]*(n^θ) - n - price*cf
        W_future = β*transpose(F[s_index, :])*W

        if W_future < 0
            W_temp[s_index] = prof
            exit[s_index] = 1
        elseif W_future >= 0
            W_temp[s_index] = prof + W_future
            exit[s_index] = 0
        end
    end
    res.pol_func = exit

    return W_temp
end

function v_iterate(par::Params, res::Results; tol::Float64=1e-3, price::Float64 = 0.7)
    err = 100.0
    while err > tol
        W_temp = firms_problem(par, res; price)
        err = maximum(abs.(res.W .- W_temp))
        res.W = W_temp
    end
end


# Firm's problem for a version with shocks

function firms_problem_shocks(par::Params, res::Results; price::Float64 = 0.7)
    @unpack s, Ns, F, β, θ, cf, α = par
    @unpack W, pol_func, U0, V0, V1 = res

    U1 = zeros(5)
    eulercon = 0.5772156649

    for s_index = 1:Ns
        n = (1/(θ*price*s[s_index]))^(1/(θ-1))
        prof = price*s[s_index]*(n^θ) - n - price*cf
        V1[s_index] = prof
        V0[s_index] = prof + β*(transpose(U0)*F[s_index, :])
        U1[s_index] = eulercon/α + 1/α*log(exp(α*V0[s_index]) + exp(α*V1[s_index]))
        res.labor[s_index] = n
    end
    return U1
end

function v_iterate_shocks(par::Params, res::Results; tol::Float64 = 1e-3, price::Float64 = 0.7)
    err = 100.0
    it = 1
    while err > tol
        U_temp = firms_problem_shocks(par, res; price)
        err = maximum(abs.(res.U0 .- U_temp))
        res.U0 = U_temp
        it = 2
    end
end

function compute_sigma(par::Params, res::Results)
    @unpack α = par
    @unpack V0, V1, σ = res

    σ = exp.(α*V1)./(exp.(α*V1) .+ exp.(α*V0))

    return σ
end


# Function to compute value of entrants
function value_entrants(par::Params, res::Results; price::Float64 = 0.7)
    @unpack α = par

    if α == 0
        value_entr = transpose(res.W)*par.ν
    else
        value_entr = transpose(res.U0)*par.ν
    end

    return value_entr
end

# Function to find an equilibrium price
function find_price(par::Params, res::Results; tol = 10e-4)
    @unpack ce, α = par

    # Guess price
    price = 0.7

    # Difference for the loop over the free entry condition
    dif = 10

    while (dif > tol)
        if α == 0
            v_iterate(par, res; price = price)
        else
            v_iterate_shocks(par, res; price = price)
        end

        val_entr = value_entrants(par, res; price = price)
        costs = price*ce
        dif = abs(val_entr - costs)

        if (dif > tol)
            if α == 1
                step = 0.0001
            elseif α == 0
                step = 0.001
            elseif α == 2
                step = 0.00001
            end
            if val_entr > costs
                price = price - step
            else
                price = price + step
            end
        #    println("updated price = $price and dif = $dif")
        end
        if (dif <= tol)
        #    println("Equilibrium price = $price")
        end
    end
    res.p = price
end


function solve_model(; cf::Int64, α::Int64)
    par, res = Initialize(; cf, α)

    @unpack α = par

    res.p = find_price(par, res)

    if α > 0
        res.pol_func = compute_sigma(par, res)
    end

    return par, res
end
