#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 3
    AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

# Finds the optimal labor supply
function getLabor(a_index, z_index, age, ap_index, par::parameters, res::results)
    @unpack R, a_grid, δ = par
    @unpack γ, θ, e, w, r = res
    
    l = (γ*(1-θ)*e[age, z_index]*w - (1-γ)*((1+r - δ)*a_grid[a_index] - a_grid[ap_index])) / ((1-θ)*w*e[age, z_index])
    
    if l > 1.0
        l = 1.0 
    elseif l < 0.0 
        l = 0.0 
    end 
    
    return l 
end

# Finds utility of retirees
function UtilityRetiree(c::Float64, σ::Float64, γ::Float64)
    if c > 0
        u_r = (c^((1-σ)*γ))/(1-σ) # CRRA utility for retiree
    else
        u_r = -Inf
    end

    u_r 
end

# Finds utility of workers
function UtilityWorker(c::Float64, l::Float64, σ::Float64, γ::Float64)
    if c > 0
        u_w = (((c^γ)*(1-l)^(1-γ))^(1-σ))/(1-σ) # CRRA utility for worker
    else
        u_w = -Inf
    end

    u_w 
end

# Bellman function for retirees
function bellman_retiree(par::parameters, res::results, age::Int64)
    @unpack a_grid, na, nz, N, R, σ, β, π = par
    @unpack val_func, pol_func, θ, γ, e, z, r, b, w = res

    for z_index = 1:nz
        for age = N:-1:R
            if age == N
                for a_index in 1:na
                    c = (1 + r) * a_grid[a_index] + b
                    pol_func[a_index, z_index, N] = 0.0

                    if c > 0
                        val_func[a_index, z_index, N] = UtilityRetiree(c, σ, γ)
                    end
                end
            else 
                choice_lower = 1                                # for exploiting monotonicity of policy function

                for a_index in 1:na
                    a = a_grid[a_index]
                    max_val = -Inf    

                    for ap_index in choice_lower:na
                        a_prime = a_grid[ap_index]
                        c = (1 + r) * a_grid[a_index] + b - a_prime

                        if c > 0
                            val = UtilityRetiree(c, σ, γ) + β*val_func[ap_index, z_index, age+1]
                            
                            if val > max_val
                                max_val = val
                                pol_func[a_index, z_index, age] = a_prime
                                choice_lower = ap_index
                            end
                        end
                    end

                    val_func[a_index, z_index, age] = max_val
                end
            end
        end
    end
end

# Bellman function for workers
function bellman_worker(par::parameters, res::results, age::Int64)
    @unpack a_grid, na, nz, N, R, σ, β, π = par
    @unpack val_func, pol_func, labor, θ, γ, e, z, w, r = res

    for z_index = 1:nz
        for age = R-1:-1:1
            choice_lower = 1                                                # for exploiting monotonicity of policy function

            for a_index in 1:na
                a = a_grid[a_index]
                max_val = -Inf    

                for ap_index in choice_lower:na
                    a_prime = a_grid[ap_index]
                    l = getLabor(a_index, z_index, age, ap_index, par, res)
                    c = w * (1-θ) * e[age, z_index] * l + (1+r) * a_grid[a_index] - a_prime

                    if c > 0 && l >= 0 && l <= 1
                        if age == R-1
                            val = UtilityWorker(c, l, σ, γ) + β*(val_func[ap_index, z_index, age+1])
                        else
                            val = UtilityWorker(c, l, σ, γ) + β*(val_func[ap_index, 1, age+1]*π[z_index,1] + val_func[ap_index, 2, age+1]*π[z_index,2])
                        end

                        if val > max_val
                            max_val = val
                            pol_func[a_index, z_index, age] = a_prime
                            labor[a_index, z_index, age] = l
                            choice_lower = ap_index
                        end
                    end
                end

                val_func[a_index, z_index, age] = max_val
            end
        end
    end
end

# Iterate backwards
function V_iterate(par::parameters, res::results)
    @unpack N, R = par

    for age in N:-1:1
        if age >= R                  
            bellman_retiree(par, res, age)    
        else                                   
            bellman_worker(par, res, age)     
        end
    end
end

# Create population distribution by age
function AgeDistribution(N::Int64, n::Float64)
    μ = ones(N)    
    μ[1] = 1      

    for i in 2:N   
        μ[i] = μ[i-1]/(1+n)
    end
    μ = μ ./sum(μ) 
    μ 
end

# Generate stationary distribution
function get_distr(par::parameters, res::results)
    @unpack π, π0, na, nz, a_grid, N, n = par 
    @unpack pol_func = res 
    
    F = zeros(par.na, par.nz, par.N)
    F[1,:, 1] = par.μ_1 .* π0  # we know what they should start with at age = 1
                                            # everyone has a = 0 at age = 1.
    for age = 2:N
        for z_index = 1:2
            kindexes = findall(!iszero,F[:,z_index, age-1])
            for k_index in kindexes
                kp = pol_func[k_index, z_index, age-1]
                kp_index = searchsortedfirst(a_grid, kp)
                F[kp_index,1, age] += F[k_index, z_index, age-1] * π[z_index, 1] / (1+n)
                F[kp_index,2, age] += F[k_index, z_index, age-1] * π[z_index, 2] / (1+n)
            end
        end
    end
    res.F = F

    return F
end

# Calculate wage, interest rate, and pension benefit using aggregate capital and labor
function CalcPrices(par::parameters, res::results, K::Float64, L::Float64)
    @unpack α, δ, R, N, n = par
    @unpack F, θ = res

    # Calculate w
    res.w = (1-α) * K^α * L^(-α)

    # Calculate r
    res.r = α * K^(α-1) * L^(1-α) - δ

    # Calculate b
    if θ == 0.0
        res.b = 0.0
    else
        @unpack w = res

        μ = AgeDistribution(N, n)
        res.b = (θ*w*L)/ sum(μ[R:N])
    end
end

# Calculate aggregate capital and labor
function CalcAggregate(par::parameters, res::results)
    @unpack R, N, na, nz, a_grid = par
    @unpack F, labor, e = res

    K1 = 0.0
    L1 = 0.0

    # Calculate aggregate capital demand (K1)
    K1 = sum(F .* a_grid)

    # Calculate aggregate labor demand (L1)
    for age = 1:R-1
        for z_index = 1:nz
            for k_index = 1:na
                L1 += (F[k_index, z_index, age] .* e[age, z_index] .* labor[k_index, z_index, age])
            end
        end
    end

    return K1, L1
end

# Calculate welfare objects - need to sum over result in order to calculate total welfare
function CalcWelfare(res::results)
    @unpack val_func, F = res

    # Calculate welfare
    welfare = val_func .* F
    welfare = sum(welfare[isfinite.(welfare)])

    return welfare
end

# Calculate CV
function computecv(par::parameters, res::results)
    @unpack a_grid, na, nz, N, R = par
    @unpack θ, w, r, b, e, labor, F = res

    mean_wealth = 0.0
    var_wealth  = 0.0
    cv_wealth   = 0.0

    for i=1:par.na, j=1:par.nz, age=1:par.N
        if age < par.R
            mean_wealth += (res.w * (1-res.θ) * e[age, j] * res.labor[i,j,age] + (1+res.r) * par.a_grid[i]) * res.F[i,j,age]
        elseif age >= par.R
            mean_wealth += ((1 + res.r) * par.a_grid[i] + res.b) * res.F[i,j,age]
        end
    end

    for i=1:par.na, j=1:par.nz, age=1:par.N
        if age < par.R
            var_wealth += (res.w * (1 - res.θ) * e[age, j] * res.labor[i,j,age] + (1+res.r) * par.a_grid[i] - mean_wealth)^2 * res.F[i,j,age]
        elseif age >= par.R
            var_wealth += ((1 + res.r) * par.a_grid[i] + res.b - mean_wealth)^2 * res.F[i,j,age]
        end
    end

    cv_wealth = var_wealth^(0.5)/mean_wealth
    return cv_wealth
end

# Solve the model
function SolveModel(;θ::Float64 = 0.11, z::Vector{Float64} = [3.0, 0.5], γ::Float64 = 0.42, λ = 0.5, iter = 1000, tol = 0.005)
    par, res = Initialize(θ, z, γ)

    # Initial guesses
    if λ == 0.3
        K0 = 1.01
        L0 = 0.17
    else
        K0 = 3.3
        L0 = 0.3
    end 
    
    diff = 10.0
    n = 1

    # Find K, L
    while (diff > tol && n < iter)
        println("BEGINNING ITERATION $n")
        n = n+1
        
        # Calculate w, r, b based on guess of K0, L0
        CalcPrices(par, res, K0, L0)

        # Solve
        V_iterate(par, res)
        get_distr(par, res)

        K1, L1 = CalcAggregate(par, res)
        diff = abs(K1 - K0) + abs(L1 - L0)
        println("FINDS DIFFERENCE $diff") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")
        println("\tK0 = $K0 \t K1 = $K1")
        println("\tL0 = $L0 \t L1 = $L1\n")

        # Adjust guess
        if diff > tol
            K0 = λ*K1 + (1-λ)*K0
            L0 = λ*L1 + (1-λ)*L0
        else
            println("DONE!\n")
            break
        end
    end

    # NEED TO OUTPUT: K, L, w, r, b, W, cv
    res.K = K0
    res.L = L0
    welfare = CalcWelfare(res)
    cv = computecv(par, res)

    # Produce results matrix
    output = [res.K; res.L; res.w; res.r; res.b; welfare; cv]

    if (θ == 0.11) && sum(z) >2.0 && (γ != 1.0 )
        return output, par, res
    else
        return output
    end
end