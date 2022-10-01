#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 3
    AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

# Finds the optimal labor supply
function getLabor(a_index, z_index, age, ap_index, prim::primitives, ins::inputs, w, r)
    @unpack e, γ, R, δ = prim
    @unpack θ,a_grid = ins 
    
    l = (γ*(1-θ)*e[age, z_index]*w - (1-γ)*((1+r - δ)*a_grid[a_index] - a_grid[ap_index])) / ((1-θ)*w*e[age, z_index])
    
    if l > 1.0
        l = 1.0 
    elseif l < 0.0 
        l = 0.0 
    end 
    
    return l 
end

# Finds utility of retirees
function UtilityRetiree(c::Float64, σ::Int64, γ::Float64)
    if c > 0
        u_r = (c^((1-σ)*γ))/(1-σ) # CRRA utility for retiree
    else
        u_r = -Inf
    end

    u_r 
end

# Finds utility of workers
function UtilityWorker(c::Float64, l::Float64, σ::Int64, γ::Float64)
    if c > 0
        u_w = (((c^γ)*(1-l)^(1-γ))^(1-σ))/(1-σ) # CRRA utility for worker
    else
        u_w = -Inf
    end

    u_w 
end

# Bellman function for retirees
function bellman_retiree(prim::primitives, ins::inputs, w, r, b, θ, val_func, pol_func)
    @unpack N, R, σ, β, π, γ, e, z = prim
    @unpack a_grid, na = ins

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
    return val_func, pol_func
end

# Bellman function for workers
function bellman_worker(prim::primitives, ins::inputs, w::Float64, r::Float64, θ::Float64)
    @unpack N, R, σ, β, π, γ, e, z, = prim
    @unpack a_grid, na = ins

    # initialize grids for val, pol, and labor
    val_func = zeros(na, 2, N)
    pol_func = zeros(na, 2, N)
    labor    = zeros(na, 2, N)

    for z_index = 1:2
        for age = R-1:-1:1
            choice_lower = 1                                                # for exploiting monotonicity of policy function

            for a_index in 1:na
                a = a_grid[a_index]
                max_val = -Inf    

                for ap_index in choice_lower:na
                    a_prime = a_grid[ap_index]
                    l = getLabor(a_index, z_index, age, ap_index, prim, ins, w, r)
                    
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
    return pol_func, val_func, labor
end

# Iterate backwards
# function V_iterate(par::parameters, res::results)
#     @unpack N, R = par

#     for age in N:-1:1
#         if age >= R                  
#             bellman_retiree(par, res, age)    
#         else                                   
#             bellman_worker(par, res, age)     
#         end
#     end
# end

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
function get_ss_distr(pol_func, prim::primitives, ins::inputs)
    @unpack π, π0, N, n = prim
    @unpack na, a_grid = ins
    
    F = zeros(na, 2, N)
    F[1,:, 1] = prim.μ_1 .* π0  # we know what they should start with at age = 1
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

    return F
end

# Calculate wage, interest rate, and pension benefit using aggregate capital and labor
# function CalcPrices(par::parameters, res::results, K::Float64, L::Float64)
#     @unpack α, δ, R, N, n = par
#     @unpack F, θ = res

#     # Calculate w
#     res.w = (1-α) * K^α * L^(-α)

#     # Calculate r
#     res.r = α * K^(α-1) * L^(1-α) - δ

#     # Calculate b
#     if θ == 0.0
#         res.b = 0.0
#     else
#         @unpack w = res

#         μ = AgeDistribution(N, n)
#         res.b = (θ*w*L)/ sum(μ[R:N])
#     end
# end

# # Calculate aggregate capital and labor
# function CalcAggregate(par::parameters, res::results)
#     @unpack R, N, na, nz, a_grid = par
#     @unpack F, labor, e = res

#     K1 = 0.0
#     L1 = 0.0

#     # Calculate aggregate capital demand (K1)
#     K1 = sum(F .* a_grid)

#     # Calculate aggregate labor demand (L1)
#     for age = 1:R-1
#         for z_index = 1:nz
#             for k_index = 1:na
#                 L1 += (F[k_index, z_index, age] .* e[age, z_index] .* labor[k_index, z_index, age])
#             end
#         end
#     end

#     return K1, L1
# end


# # Solve the model
# function SolveModel(;θ::Float64 = 0.11, z::Vector{Float64} = [3.0, 0.5], γ::Float64 = 0.42, λ = 0.5, iter = 1000, tol = 0.005)
#     par, res = Initialize(θ, z, γ)

#     # Initial guesses
#     if λ == 0.3
#         K0 = 1.01
#         L0 = 0.17
#     else
#         K0 = 3.3
#         L0 = 0.3
#     end 
    
#     diff = 10.0
#     n = 1

#     # Find K, L
#     while (diff > tol && n < iter)
#         println("BEGINNING ITERATION $n")
#         n = n+1
        
#         # Calculate w, r, b based on guess of K0, L0
#         CalcPrices(par, res, K0, L0)

#         # Solve
#         V_iterate(par, res)
#         get_distr(par, res)

#         K1, L1 = CalcAggregate(par, res)
#         diff = abs(K1 - K0) + abs(L1 - L0)
#         println("FINDS DIFFERENCE $diff") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")
#         println("\tK0 = $K0 \t K1 = $K1")
#         println("\tL0 = $L0 \t L1 = $L1\n")

#         # Adjust guess
#         if diff > tol
#             K0 = λ*K1 + (1-λ)*K0
#             L0 = λ*L1 + (1-λ)*L0
#         else
#             println("DONE!\n")
#             break
#         end
#     end

#     # NEED TO OUTPUT: K, L, w, r, b, W, cv
#     res.K = K0
#     res.L = L0

# end