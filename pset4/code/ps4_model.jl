### more functions here ###

function initialGuesses(ins::inputs, res::results)
    @unpack  T = ins
    @unpack  K, L = res

    # Guess for K and L - T points from K[1] (L[1]) and K[end] (L[end]) equally spaced
    K_step = (K[end] - K[1])/(T-1)
    L_step = (L[end] - L[1])/(T-1)
    for i = 2:T-1
        K[i] = K[i-1] + K_step
        L[i] = L[i-1] + L_step
    end

    # Path of wages w
#    res.w = (1-α) * K.^α .* L.^(-α)

    # Path of interest rates r
    #res.r = α * K.^(α-1) .* L.^(1-α) .- δ

   # Path of pensions
 #  for t = 1:T
#        if t < t_noss
#            @unpack w = res
#            μ = AgeDistribution(N, n)
#            res.b[t] = (θ*w[t]*L[t])/ sum(μ[R:N])
#        else
#            res.b[t] = 0
#        end
#    end
    return K, L
end

function prices(prim::primitives, ins::inputs, res::results; K0::Array{Float64}, L0::Array{Float64})
    @unpack α, δ, R, N, n = prim
    @unpack θ, T, t_noss = ins
    @unpack K, L, Γ = res

    # Path of wages w
    res.w = (1-α) * K.^α .* L.^(-α)

    # Path of interest rates r
    res.r = α * K.^(α-1) .* L.^(1-α) .- δ

    # Path of pensions
    for t = 1:T
        if t < t_noss
            @unpack w = res
            μ = AgeDistribution(N, n)
            res.b[t] = (θ*w[t]*L[t])/ sum(μ[R:N])
        else
            res.b[t] = 0
        end
    end
    return res.w, res.r, res.b
end

# Bellman function for retirees
function bellman_retiree_time(prim::primitives, ins::inputs,  val_func, pol_func)
    @unpack N, R, σ, β, π, γ, e, z = prim
    @unpack a_grid, na, T, t_noss, θ  = ins
    @unpack w, r, b = res


    for t = T:-1:1
    #    println("Period t = $t")

        wage = w[t]
        rate = r[t]
        ben  = b[t]
        for z_index = 1:2
            for age = N:-1:R
        #        println("Age  = $age")

                if age == N
                    for a_index in 1:na
                        c = (1 + rate) * a_grid[a_index] + ben
                        pol_func[t, a_index, z_index, N] = 0.0

                        if c > 0
                            val_func[t, a_index, z_index, N] = UtilityRetiree(c, σ, γ)
                        end
                    end
                else
                    choice_lower = 1                                # for exploiting monotonicity of policy function

                    for a_index in 1:na
                        a = a_grid[a_index]
                        max_val = -Inf

                        for ap_index in choice_lower:na
                            a_prime = a_grid[ap_index]
                            c = (1 + rate) * a_grid[a_index] + ben - a_prime

                            if c > 0
                                if t < T
                                    val = UtilityRetiree(c, σ, γ) + β*val_func[t+1, ap_index, z_index, age+1]
                                else
                                    val = val_func[end,  ap_index, z_index, age+1]
                                end
                                if val > max_val
                                    max_val = val
                                    pol_func[t, a_index, z_index, age] = a_prime
                                    choice_lower = ap_index
                                end
                            end
                        end

                        val_func[t, a_index, z_index, age] = max_val
                    end
                end
            end
        end
    end

    return val_func, pol_func
end

# Bellman function for workers
function bellman_worker_time(prim::primitives, ins::inputs, res::results)
    @unpack N, R, σ, β, π, γ, e, z, = prim
    @unpack a_grid, na, T, t_noss, θ = ins
    @unpack w, r, b, val_func, pol_func, labor = res


    for t = T:-1:1
    #    println("Period t = $t")
        wage = w[t]
        rate = r[t]
        ben  = b[t]
        for z_index = 1:2
            for age = R-1:-1:1
            #    println("Age  = $age")

                choice_lower = 1                                                # for exploiting monotonicity of policy function

                for a_index in 1:na
                    a = a_grid[a_index]
                    max_val = -Inf

                    for ap_index in choice_lower:na
                        a_prime = a_grid[ap_index]
                        l = getLabor(a_index, z_index, age, ap_index, prim, ins, wage, rate)

                        c = wage * (1-θ) * e[age, z_index] * l + (1+rate) * a_grid[a_index] - a_prime

                        if c > 0 && l >= 0 && l <= 1
                            if t < T
                                if age == R-1
                                    val = UtilityWorker(c, l, σ, γ) + β*(val_func[t+1, ap_index, z_index, age+1])
                                else
                                    val = UtilityWorker(c, l, σ, γ) + β*(val_func[t+1, ap_index, 1, age+1]*π[z_index,1] + val_func[t+1, ap_index, 2, age+1]*π[z_index,2])
                                end
                            else
                                val = val_func[end, ap_index, z_index, age+1]
                            end

                            if val > max_val
                                max_val = val
                                pol_func[t, a_index, z_index, age] = a_prime
                                labor[t, a_index, z_index, age] = l
                                choice_lower = ap_index
                            end
                        end
                    end

                    val_func[t, a_index, z_index, age] = max_val
                end
            end
        end
    end
    return val_func, pol_func, labor
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
function get_ss_distr_time(pol_func, prim::primitives, ins::inputs)
    @unpack π, π0, N, n = prim
    @unpack na, a_grid, T = ins


    F = zeros(T, na, 2, N)
    for t = 1:T
        F[t, 1,:, 1] = prim.μ_1 .* π0  # we know what they should start with at age = 1
    end
    for t = 1:T
        for age = 2:N
            for z_index = 1:2
                kindexes = findall(!iszero,F[t, :,z_index, age-1])
                for k_index in kindexes
                    kp = pol_func[t, k_index, z_index, age-1]
                #    kp_index = searchsortedfirst(a_grid, kp)
                    kp_index = findall(x -> x == kp, a_grid)[1]
                    F[t, kp_index,1, age] += F[t, k_index, z_index, age-1] * π[z_index, 1] / (1+n)
                    F[t, kp_index,2, age] += F[t, k_index, z_index, age-1] * π[z_index, 2] / (1+n)
                end
            end
        end
    end
    return F
end

function CalcAggregate_time(prim::primitives, ins::inputs, res::results)
     @unpack R, N, e = prim
     @unpack Γ, labor = res
     @unpack T, na, a_grid = ins


     K1 = zeros(T)
     L1 = zeros(T)

     # Calculate aggregate capital demand (K1)

     # Calculate aggregate labor demand (L1)
    for t = 1:T
        K1[t] = sum(Γ[t, :, :, :] .* ins.a_grid)
        for age = 1:R-1
            for z_index = 1:2
                for k_index = 1:na
                    L1[t] += (Γ[t, k_index, z_index, age] .* e[age, z_index] .* labor[t, k_index, z_index, age])
                end
            end
        end
    end

    return K1, L1
end

### solveModel
function solveModel(prim::primitives, ins::inputs, res::results; tol = 0.001, iter = 1000, λ = 0.95)
    @unpack N, R, σ, β, π, γ, e, z, = prim
    @unpack a_grid, na, T, t_noss, θ = ins
    @unpack w, r, b, labor, pol_func, val_func = res

    # Initial guesses for capital and labor paths along the transition
    K0, L0 = initialGuesses(ins, res)

    dif = 100
    n_iter = 1

    # Solve the agent's problem
    while (dif > tol && n_iter < iter)

        # Prices paths
        w, r, b = prices(prim, ins, res; K0 = K0, L0 = L0)

        val_worker, pol_worker, res.labor = bellman_worker_time(prim, ins, res)
        val_func, pol_func = bellman_retiree_time(prim, ins,  val_worker, pol_worker)

        # Find stationary distribution
        res.Γ = get_ss_distr_time(pol_func,  prim, ins)
        # Find aggregate variables according to the obtained decision rule and distribution
        K1, L1 = CalcAggregate_time(prim, ins, res)

        # Compare
        max_dif_K = maximum(abs.(K1 .- K0))
        max_dif_L = maximum(abs.(L1 .- L0))

        dif = maximum([max_dif_K, max_dif_L])
        println("FINDS DIFFERENCE $dif") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")
        println("\tK0 = $K0 \t K1 = $K1")
        println("\tL0 = $L0 \t L1 = $L1\n")

        if dif > tol
            K0 = λ*K0 + (1-λ)*K1
            L0 = λ*L0 + (1-λ)*L1
        else
            println("DONE!\n")
            break
        end

        n_iter += 1
    end

    ##############

    return ins, res     # what we need in order to make graphs
end
