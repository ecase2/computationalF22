#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 3
    AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

# construct the grid for consumption and labor
# it should make the code faster
function fill_end_grids(par::parameters, res::results)
    @unpack a_grid, na, nz, N, R = par
    @unpack w, r, b, z, θ, γ, e = res

    for a_index = 1:na, z_index = 1:nz, j = 1:N, ap_index = 1:na
        if j < R
            res.l_grid[a_index, z_index, j, ap_index] = (γ*(1-θ)*e[j, z_index]*w - (1-γ)*((1+r)*a_grid[a_index] - a_grid[ap_index])) / ((1-θ)*w*e[j, z_index])
            if res.l_grid[a_index, z_index, j, ap_index] < 0
                res.l_grid[a_index, z_index, j, ap_index] = 0
            end
            if res.l_grid[a_index, z_index, j, ap_index] > 1
                res.l_grid[a_index, z_index, j, ap_index] = 1
            end
            res.c_grid[a_index, z_index, j, ap_index] = w*(1-θ)*e[j, z_index]*res.l_grid[a_index, z_index, j, ap_index] + (1+r)*a_grid[a_index] - a_grid[ap_index]
            # important: le may be negative according to the formula - we will need to proceed only if l is between 0 and 1.
        end
        if j >= R
            #grid.l[a_index, z_index, j, ap_index] = 0 # redundant
            res.c_grid[a_index, z_index, j, ap_index] = (1+r)*a_grid[a_index] + b - a_grid[ap_index]
        end
    end
end

# Do backwards induction to get value function and policy function
function backward_iteration(par::parameters, res::results)
    @unpack a_grid, na, nz, N, R, σ, β, π = par
    @unpack val_func, pol_func, c_grid, l_grid, θ, γ, e, z = res # w, r, b, 

    # last period - save nothing, consume everythinig
    for a_index = 1:na, z_index = 1:nz
        pol_func[a_index, z_index, N] = 0
        cons                          = c_grid[a_index, z_index, N, 1]
        val_func[a_index, z_index, N] = cons^((1-σ)*γ)/(1-σ)
    end

    t = N-1
    while t > 0
        # println("age is ", t)
        if t >= R
            for a_index = 1:na, z_index = 1:nz # not very efficient. Make to do nz loops which are identical
                max_val = -Inf
                for ap_index = 1:na
                    cons = c_grid[a_index, z_index, t, ap_index]
                    if cons >= 0
                        util = cons^((1-σ)*γ)/(1-σ) + β*val_func[ap_index, z_index, t + 1]
                        if util > max_val
                            max_val = util
                            res.pol_func[a_index, z_index, t] = a_grid[ap_index]
                            res.val_func[a_index, z_index, t] = util
                        end
                    end
                end
            end
        end
        if t < R
            for a_index = 1:na, z_index = 1:nz
                max_val = -Inf
                for ap_index = 1:na
                    lab = l_grid[a_index, z_index, t, ap_index]  # CHANGED FROM L TO LGRID
                    cons = c_grid[a_index, z_index, t, ap_index] # CHANGED FROM C TO CGRID
                    if cons >= 0
                        util = (cons^γ*(1-lab)^(1-γ))^(1-σ)/(1-σ) + β*(transpose(π[z_index, :])*val_func[ap_index, :, t + 1])
                        if util > max_val
                            max_val = util
                            res.pol_func[a_index, z_index, t] = a_grid[ap_index]
                            res.val_func[a_index, z_index, t] = util
                            res.l[a_index, z_index, t] = lab
                        end
                    end
                end
            end

        end
        t -=  1
    end
end

# Generate stationary distribution
function get_distr(par::parameters, res::results)
    @unpack π, π0, na, a_grid, N, n = par 
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
    @unpack α, δ, R, N = par
    @unpack F, θ = res

    # Calculate w
    res.w = (1-α)*(K^α)*(L^(-α))

    # Calculate r
    res.r = α*(K^(α-1))*L^(1-α) - δ

    # Calculate b
    @unpack w = res
    res.b = (θ*w*L)/ sum(F[ :, :, R:N])
end

# Determine market clearing by calculating aggregate capital demand and aggregate labor demand
function AggregateDemand(par::parameters, res::results)
    @unpack R, N, na, nz = par
    @unpack F, l, pol_func, e = res

    K1 = 0.0
    L1 = 0.0

    # Calculate aggregate capital demand (K1)
    for age = 1:N
        for z_index = 1:nz
            for k_index = 1:na
                K1 += (F[k_index, z_index, age] .* pol_func[k_index, z_index, age])
            end
        end
    end

    # Calculate aggregate labor demand (L1)
    for age = 1:R-1
        for z_index = 1:nz
            for k_index = 1:na
                L1 += (F[k_index, z_index, age] .* e[age, z_index] .* l[k_index, z_index, age])
            end
        end
    end

    return K1, L1
end

# Calculate welfare
function CalcWelfare(res::results)
    @unpack val_func, F = res
    
    # Calculate welfare
    welfare = val_func .* F
    welfare = sum(welfare)

    return welfare
end

# Calculate CV
function CalcCV(res::results)
    @unpack val_func, F = res

    # Calculate welfare
    welfare = CalcWelfare(res)

    # Calculate mean welfare
    welfare_mean = mean(welfare)

    # Calculate standard deviation
    welfare_sd = std(welfare)

    # Calculate CV
    cv = welfare_sd/welfare_mean

    return cv
end

# Solve the model
function SolveModel(;θ::Float64 = 0.11, z::Vector{Float64} = [3.0, 0.5], γ::Float64 = 0.42, iter = 1000, tol = 0.005)
    par, res = Initialize(θ, z, γ)
   
    # Set initial guesses
    K0 = 3.3
    L0 = 0.3
    diff = 10.0
    n = 1

    # Find K, L
    while (diff > tol && n < iter)
        println("BEGINNING ITERATION $n")
        n = n+1
        
        # Solve
        fill_end_grids(par, res)
        backward_iteration(par, res)
        get_distr(par, res)

        # Calculate w, r, b based on guess of K0, L0
        CalcPrices(par, res, K0, L0)

        K1, L1 = AggregateDemand(par, res) 
        diff = abs(K1 - K0) + abs(L1 - L0) 
        println("FINDS DIFFERENCE $diff... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")


        # Adjust guess
        if diff > tol 
            if (K1 - K0) > 0 # if the demand is higher than supply, new supply guess has more weight
                λ_k = 0.3
            else 
                λ_k = 0.7
            end 
            if (L1 - L0) >0 
                λ_l = 0.3 
            else 
                λ_l = 0.7
            end 
            K0 = λ_k*K1 + (1-λ_k)*K0 
            L0 = λ_l*L1 + (1-λ_l)*L0 
        else 
            break
        end 
    end

    
    # NEED TO OUTPUT: K, L, w, r, b, W, cv 
    res.K = K0
    res.L = L0
    welfare = CalcWelfare(res)
    cv = CalcCV(res)

    # Produce results matrix
    output = [res.K; res.L; res.w; res.r; res.b; welfare; cv]
    
    return output    
end