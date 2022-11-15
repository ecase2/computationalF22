#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 3
    AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

# Finds the optimal labor supply
function getLabor(a_index, z_index, age, ap_index, prim::primitives, res::StationaryResults)
    @unpack R, a_grid, δ, γ, θ, e = prim
    @unpack w, r = res

    l = (γ*(1-θ)*e[age, z_index]*w - (1-γ)*((1+r)*a_grid[a_index] - a_grid[ap_index])) / ((1-θ)*w*e[age, z_index])

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
function bellman_retiree(prim::primitives, res::StationaryResults, age::Int64)
    @unpack a_grid, na,  N, R, σ, β, π, θ, γ, e, z= prim
    @unpack val_func, pol_func, r, b, w = res

    for z_index = 1:2
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
function bellman_worker(prim::primitives, res::StationaryResults, age::Int64)
    @unpack a_grid, na,  N, R, σ, β, π, θ, γ, e, z = prim
    @unpack val_func, pol_func, labor, w, r = res

    for z_index = 1:2
        for age = R-1:-1:1
            choice_lower = 1                                                # for exploiting monotonicity of policy function

            for a_index in 1:na
                a = a_grid[a_index]
                max_val = -Inf

                for ap_index in choice_lower:na
                    a_prime = a_grid[ap_index]
                    l = getLabor(a_index, z_index, age, ap_index, prim, res)
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
function V_iterate(prim::primitives, res::StationaryResults)
    @unpack N, R = prim

    for age in N:-1:1
        if age >= R
            bellman_retiree(prim, res, age)
        else
            bellman_worker(prim, res, age)
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
function get_distr(prim::primitives, res::StationaryResults)
    @unpack π, π0, na, a_grid, N, n = prim
    @unpack pol_func = res

    F = zeros(na, 2, N)
    F[1,:, 1] = prim.μ_1 .* π0  # we know what they should start with at age = 1
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
    res.Γ = F

    return res.Γ
end

# Calculate wage, interest rate, and pension benefit using aggregate capital and labor
function CalcPrices(prim::primitives, res::StationaryResults, K::Float64, L::Float64)
    @unpack α, δ, R, N, n,  θ  = prim
    @unpack Γ = res

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
function CalcAggregate(prim::primitives, res::StationaryResults)
    @unpack R, N, na, a_grid, e = prim
    @unpack Γ, labor = res

    K1 = 0.0
    L1 = 0.0

    # Calculate aggregate capital demand (K1)
    K1 = sum(Γ .* a_grid)

    # Calculate aggregate labor demand (L1)
    for age = 1:R-1
        for z_index = 1:2
            for k_index = 1:na
                L1 += (Γ[k_index, z_index, age] .* e[age, z_index] .* labor[k_index, z_index, age])
            end
        end
    end

    return K1, L1
end

# Calculate welfare objects - need to sum over result in order to calculate total welfare
function CalcWelfare(res::StationaryResults)
    @unpack val_func, Γ = res

    # Calculate welfare
    welfare = val_func .* Γ
    welfare = sum(welfare[isfinite.(welfare)])

    return welfare
end

# Calculate CV
function computecv(prim::primitives, res::StationaryResults)
    @unpack θ, e, a_grid, na,  N, R = prim
    @unpack w, r, b, labor, Γ = res

    mean_wealth = 0.0
    var_wealth  = 0.0
    cv_wealth   = 0.0

    for i=1:prim.na, j=1:2, age=1:prim.N
        if age < prim.R
            mean_wealth += (res.w * (1-prim.θ) * e[age, j] * res.labor[i,j,age] + (1+res.r) * prim.a_grid[i]) * res.Γ[i,j,age]
        elseif age >= prim.R
            mean_wealth += ((1 + res.r) * prim.a_grid[i] + res.b) * res.Γ[i,j,age]
        end
    end

    for i=1:prim.na, j=1:2, age=1:prim.N
        if age < prim.R
            var_wealth += (res.w * (1 - prim.θ) * e[age, j] * res.labor[i,j,age] + (1+res.r) * prim.a_grid[i] - mean_wealth)^2 * res.Γ[i,j,age]
        elseif age >= prim.R
            var_wealth += ((1 + res.r) * prim.a_grid[i] + res.b - mean_wealth)^2 * res.Γ[i,j,age]
        end
    end

    cv_wealth = var_wealth^(0.5)/mean_wealth
    return cv_wealth
end

# Solve the model
function SolveModel(; al = 0.0, au = 20.0, na = 200, θ = 0.11, z = [3.0, 0.5], γ = 0.42, T = 1, t_noss = 2, λ = 0.5, iter = 1000, tol = 0.005)

    prim, res = Initialize(al, au, na, θ, z, γ, T, t_noss)

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
        CalcPrices(prim, res, K0, L0)

        # Solve
        V_iterate(prim, res)
        get_distr(prim, res)

        K1, L1 = CalcAggregate(prim, res)
        diff = abs(K1 - K0) + abs(L1 - L0)
        println("FINDS DIFFERENCE $diff") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")
        println("\tK0 = $K0 \t K1 = $K1")
        println("\tL0 = $L0 \t L1 = $L1\n")

        # Adjust guess
        if diff > tol
            K0 = λ*K1 + (1-λ)*K0
            L0 = λ*L1 + (1-λ)*L0
        else
            res.K = K1
            res.L = L1
            println("DONE!\n")
            break
        end
    end

    return  res
end
