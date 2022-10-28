#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 6
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

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
    γ_E::Float64        = 0.5772156649                           # Euler constant
    m_lb::Float64       = 1e-4                                   # lower bound on mass of entrants
    m_ub::Float64       = 10.0                                   # upper bound on mass of entrants
end

mutable struct Results
    p::Float64                                                   # price level
    m_entr::Float64                                              # mass of entrants
    L::Float64                                                   # aggregate labor
    μ::Array{Float64}                                            # distribution of firms
    pol_func::Array{Float64}                                     # exit policy function
    W::Array{Float64}                                            # value function for baseline case
    labor::Array{Float64}                                        # optimal labor demand
end

# Function to initialize parameters and results structures
function Initialize(; cf::Int64, α::Int64)
    par = Params(cf = cf, α = α)

    p = 0.5
    m_entr = 2.75
    L = 1

    μ        = ones(par.Ns)/par.Ns
    pol_func = zeros(par.Ns)
    W        = zeros(par.Ns)
    labor    = zeros(par.Ns)

    res = Results(p, m_entr, L, μ, pol_func, W, labor)

    return par, res
end

# Firm's problem for a baseline version
function firms_problem(par::Params, res::Results; price::Float64)
    @unpack s, Ns, F, β, θ, cf = par
    @unpack W, pol_func, labor = res

    W_temp = zeros(Ns)
    exit = zeros(Ns)

    for s_index = 1:Ns
        n = max((1/(θ*price*s[s_index]))^(1/(θ-1)),0)          # solve firm's static labor problem
        labor[s_index] = n

        profit = price*s[s_index]*(n^θ) - n - price*cf    # calculate profit
        W_future = β*transpose(F[s_index, :])*W         # calculate firm's continuation value

        if W_future < 0
            W_temp[s_index] = profit
            exit[s_index] = 1
        elseif W_future >= 0
            W_temp[s_index] = profit + W_future
            exit[s_index] = 0
        end
    end
    res.pol_func = exit

    return W_temp
end

# Firm's problem for a version with shocks
function firms_problem_shocks(par::Params, res::Results; price::Float64 = 0.7)
    @unpack s, Ns, F, β, θ, cf, α, γ_E = par
    @unpack W, pol_func, labor = res

    W_temp = zeros(Ns)
    exit = zeros(Ns)

    for s_index = 1:Ns
        n = max((1/(θ*price*s[s_index]))^(1/(θ-1)),0)       # solve firm's static labor problem
        labor[s_index] = n
        
        profit = price*s[s_index]*(n^θ) - n - price*cf      # calculate profit
        
        V_stay = profit + β*transpose(F[s_index, :])*W      # calculate firm's value if they stay in the market
        V_exit = profit                                     # calculate firm's value if they exit the market

        W_temp[s_index] = γ_E/α + 1/α*log(exp(α*V_stay) + exp(α*V_exit)) # calculate firm's utility
        exit[s_index] = compute_sigma(par, V_stay, V_exit)               # calculate probability of exit
    end
    res.pol_func = exit

    return W_temp
end

# Solve firm's problem
function W_iterate(par::Params, res::Results; tol::Float64=1e-3, price::Float64 = 0.7)
    @unpack α = par
    
    err = 100.0

    while err > tol
        if α == 0
            W_temp = firms_problem(par, res; price)
        else
            W_temp = firms_problem_shocks(par, res; price)
        end 

        err = maximum(abs.(res.W .- W_temp))
        res.W = W_temp
    end
end

# Compute probability that firm exits the market
function compute_sigma(par::Params, V_stay::Float64, V_exit::Float64)
    @unpack α = par

    σ = exp.(α*V_stay)./(exp.(α*V_stay) .+ exp.(α*V_exit))

    return σ
end

# Function to compute value of entrants
function value_entrants(par::Params, res::Results)
    value_entr = transpose(res.W)*par.ν

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
        W_iterate(par, res; price = price)

        val_entr = value_entrants(par, res)
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
            println("updated price = $price and dif = $dif")
        else
            println("Equilibrium price = $price")
            break
        end
    end
    res.p = price
end

# T star operator
function T_star(par::Params, res::Results)
    @unpack F, Ns, ν = par

    temp_dist = zeros(Ns)

    for s_today = 1:Ns
        for s_tomorrow = 1:Ns
            temp_dist[s_tomorrow] += (1 - res.pol_func[s_today]) * F[s_today, s_tomorrow] * res.μ[s_today]

            temp_dist[s_tomorrow] += (1 - res.pol_func[s_today]) * F[s_today, s_tomorrow] * ν[s_today] * res.m_entr
        end
    end

    return temp_dist
end

# Find stationary distribution
function find_dist(par::Params, res::Results; iter = 100, tol = 10e-4)
    @unpack pol_func, μ, m_entr = res

    diff = 10
    n = 1

    while (diff > tol && n < iter)
        temp_dist = T_star(par, res)

        diff = maximum(abs.(temp_dist .- μ))
        
        res.μ = temp_dist
        n = n + 1       
    end
end

# Calculate aggregate labor demand
function labor_demand(par::Params, res::Results)
    @unpack labor, μ, m_entr = res

    L = sum(labor .* μ) + m_entr*sum(labor .* par.ν)

    return L
end

# Calculate aggregate labor supply
function labor_supply(par::Params, res::Results)
    @unpack A, Ns, s, θ, cf, ν = par
    @unpack p, labor, μ, m_entr = res

    profit = zeros(Ns)

    for s_index = 1:Ns
        profit[s_index] = p*s[s_index]*(labor[s_index]^θ) - labor[s_index] - p*cf
    end

    agg_profit = sum(profit .* μ) + m_entr*sum(profit .* ν)

    L = 1/A - agg_profit

    return L
end

# Find the mass of entrants
function find_mass(par::Params, res::Results; iter = 100, tol = 10e-4)

    # Initialize
    m_low = par.m_lb
    m_high = par.m_ub
    res.m_entr = (m_high + m_low)/2

    # Difference for the loop over the free entry condition
    diff = 10
    n = 1

    while (diff > tol && n < iter)
        find_dist(par, res)

        L_demand = labor_demand(par, res)
        L_supply = labor_supply(par, res)

        diff = abs(L_demand - L_supply)
        
        if (diff > tol)
            if L_demand > L_supply
                m_high = res.m_entr
            else
                m_low = res.m_entr
            end

            res.m_entr = (m_high + m_low)/2
            n += 1
            println("Difference = ", diff)
        else 
            res.L = L_demand
            println("Found distribution!")
            break
        end
    end
end

# Solve model
function solve_model(; cf::Int64, α::Int64)
    par, res = Initialize(; cf, α)

    @unpack α, Ns = par

    # Run functions 
    find_price(par, res)
    find_mass(par, res)

    # Calculate output values
    mass_incumbents = transpose(1 .- res.pol_func) * res.μ
    mass_exits = transpose(res.pol_func) * res.μ

    labor_incumbents = sum(res.labor .* res.μ)
    labor_entrants = res.m_entr*sum(res.labor .* par.ν)

    res.L = labor_incumbents + labor_entrants
    frac_labor = labor_entrants / res.L

    output = [res.p; mass_incumbents; res.m_entr; mass_exits; res.L; labor_incumbents; labor_entrants; frac_labor]

    return par, res, output
end