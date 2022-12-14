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
    return K, L
end

function prices(prim::primitives, ins::inputs, res::results; Kg::Array{Float64}, Lg::Array{Float64})
    @unpack α, δ, R, N, n = prim
    @unpack θ, T, t_noss = ins
    #@unpack K, L  = res

    # Path of wages w
    for t = 1:T
        res.w[t] = (1-α) * Kg[t]^α * Lg[t]^(-α)

        # Path of interest rates r
        res.r[t] = α * Kg[t]^(α-1) * Lg[t]^(1-α) - δ
    end

    # Path of pensions
    for t = 1:T
        if t < t_noss
            @unpack w = res
            μ = AgeDistribution(N, n)
            res.b[t] = (θ*w[t]*Lg[t])/ sum(μ[R:N])
        else
            res.θ_path[t] = 0
            res.b[t] = 0
        end
    end
    return res.w, res.r, res.b
end

# Finds the optimal labor supply
function getLabor(a_index, z_index, ap_index, prim::primitives, ins::inputs, t::Int64, age::Int64)
    @unpack e, γ, R, δ = prim
    @unpack θ, a_grid = ins
    @unpack w, r, b, θ_path = res

    l = (γ*(1-θ_path[t])*e[age, z_index]*w[t] - (1-γ)*((1+r[t])*a_grid[a_index] - a_grid[ap_index])) / ((1-θ_path[t])*w[t]*e[age, z_index])

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


function optimal_labor(prim::primitives, res::results, age::Int64, a::Float64, z::Float64, ap::Float64, t::Int64) #add time
    @unpack γ, e = prim
    @unpack w, r, θ_path = res
    @unpack θ = ins


    labor_opt =  γ*(1-θ_path[t])*z*e[age]*w[t]-(1- γ)*((1+r[t])*a-ap)/((1-θ_path[t])*w[t]*z*e[age])

    #adjust labor opt if out of bounds of [0,1]
    if labor_opt > 1.0
        labor_opt = 1.0
    elseif labor_opt < 0.0
        labor_opt = 0.0
    end
    labor_opt #deliverable
end

#utility grid for retirees
function utility_grid_ret(prim::primitives, ins::inputs, res::results)
    @unpack e, γ, R, δ, σ = prim
    @unpack θ, a_grid, T, na = ins
    @unpack w, r, b, θ_path = res

    ugrid_ret = zeros(na, na, T) #one for each a and a' and t

    for t = 1:T, index_a=1:na, index_ap=1:na
        a, ap = a_grid[index_a], a_grid[index_ap]
        c = (1+r[t])*a+b[t]-ap
        if c>0
            ugrid_ret[index_a,index_ap,t] = c^((1-σ)*γ)/(1-σ)
        else
            ugrid_ret[index_a,index_ap,t] = -Inf
        end
    end
    ugrid_ret
end

#utility grid for workers
function utility_grid_work(prim::primitives, ins::inputs, res::results)
    @unpack N, R, σ, β, π, γ, e, z = prim
    @unpack a_grid, na, t_noss, θ, T  = ins
    @unpack w, r, b, θ_path = res

    ugrid_work = zeros(na*2, na, R, T) #current a and z x choice of ap x age x time
    # 4 dimension u grid, to keep track of age until retirement (R) and time

    #loop over time, current a and z, age, and  a'
    for t = 1:T, index_a = 1:na, age = 1:(R-1), index_ap = 1:na, index_z = 1:2
        a = a_grid[index_a]
        Z = z[index_z]
        ap = a_grid[index_ap]
        labor = optimal_labor(prim, res, age, a, Z, ap, t) #use optimal labor function
        c = w[t]*(1-θ_path[t])*e[age]*Z*labor+(1+r[t])*a - ap

        if c>0
            inside=(c^γ)*(1-labor)^(1-γ)
            ugrid_work[index_a+na*(index_z-1),index_ap,age,t] = (inside^(1-σ))/(1-σ)
        else
            ugrid_work[index_a+na*(index_z-1),index_ap,age,t] = -Inf
        end
    end
    ugrid_work #deliverable
end

# Bellman function for retirees - let make age an argument of the function - maybe it will make the code faster
function bellman_retiree(prim::primitives, ins::inputs, res::results, t::Int64, age::Int64, ugrid_ret::Array{Float64,3})
    @unpack N, R, σ, β, π, γ, e, z = prim
    @unpack a_grid, na, t_noss, θ  = ins
    @unpack w, r, b, pol_func, val_func = res


    value_temp = zeros(na, 2)      # value function guess

    wage = w[t]
    rate = r[t]
    ben = b[t]

        for z_index = 1:2, a_index = 1:na
            a = a_grid[a_index]
            max_val = -Inf
            if age == N
                max_val = ugrid_ret[a_index, 1, t]
            else
                for ap_index in 1:na
                    a_prime = a_grid[ap_index]
                    val = ugrid_ret[a_index, ap_index, t] + β*val_func[t+1, ap_index, z_index, age+1]
                    if val > max_val
                        max_val = val
                        pol_func[t, a_index, z_index, age] = a_prime
                    end
                end
            end
            value_temp[a_index, z_index] = max_val
        end
    return value_temp
end

# Bellman function for workers
function bellman_worker(prim::primitives, ins::inputs, res::results, t::Int64, age::Int64, ugrid_work::Array{Float64,4})
    @unpack N, R, σ, β, π, γ, e, z, = prim
    @unpack a_grid, na, T, t_noss, θ = ins
    @unpack w, r, b, pol_func, labor, val_func = res


    value_temp = zeros(na, 2)      # value function guess

    wage = w[t]
    rate = r[t]
    ben = b[t]

    for z_index = 1:2, a_index = 1:na
        a = a_grid[a_index]
        max_val = -Inf

        for ap_index in 1:na
            a_prime = a_grid[ap_index]
            val = ugrid_work[a_index+na*(z_index-1),ap_index,age,t]+β*sum(π[z_index,:].*val_func[t+1, ap_index,:,age+1])
            if val > max_val
                max_val = val
                pol_func[t, a_index, z_index, age] = a_prime
                labor[t, a_index, z_index, age] = optimal_labor(prim, res, age, a, z[z_index], a_prime, t)
            end
        end
        value_temp[a_index, z_index] = max_val
        end

    return value_temp
end

#backward induction protocol
function backward_induct(prim::primitives, ins::inputs, res::results, t::Int64, ugrid_work::Array{Float64,4}, ugrid_ret::Array{Float64,3})
    @unpack R, N = prim
    @unpack T = ins

    #loop over all ages
    for j = 1:N
        age = N-j+1 #now going backwards = 1 when j = J
        if age>= R
            res.val_func[t,:,:,age] = bellman_retiree(prim, ins, res, t, age, ugrid_ret)
        else
            res.val_func[t,:,:,age] = bellman_worker(prim, ins, res, t, age, ugrid_work)
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
function get_ss_distr(prim::primitives, ins::inputs, res::results, initial_distr::Array{Float64})
    @unpack π, π0, N, n = prim
    @unpack na, a_grid, T = ins
    @unpack pol_func = res

    distr = zeros(T, na, 2, N)
    distr[1, :, :, :] = initial_distr

    for t = 2:T
        distr[t, 1, :, 1] = prim.μ_1 .* π0  # we know what they should start with at age = 1
        for age = 2:N
            for z_index = 1:2, a_index = 1:na
                    ap = pol_func[t-1, a_index, z_index, age-1]
                    for ap_index = 1:na, zp_index = 1:2 #loop over tomorrow's states
                        if a_grid[ap_index] == ap
                            distr[t, ap_index, zp_index, age] += distr[t-1, a_index, z_index, age-1]*
                                π[z_index, zp_index]*(1/(1+n))
                        end
                    end
            end
        end
    end
    res.Γ = distr
end

function CalcAggregate(prim::primitives, ins::inputs, res::results)
     @unpack R, N, e = prim
     @unpack Γ, labor, K, L = res
     @unpack T, na, a_grid = ins


     K1 = zeros(T)
     L1 = zeros(T)

     K1[1] = K[1]
     L1[1] = L[1]

     # Calculate aggregate capital demand (K1)
     # Calculate aggregate labor demand (L1)
    for t = 2:T
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
    @unpack N, R, σ, β, π, γ, e, z = prim
    @unpack a_grid, na, T, t_noss, θ = ins
    @unpack labor, pol_func, val_func, K, L = res

    # Initial guesses for capital and labor paths along the transition
    K0, L0 = initialGuesses(ins, res)

    dif = 100
    n_iter = 1

    # Solve the agent's problem
    while (dif > tol && n_iter < iter)
        println("Iter = $n_iter")

        # Prices paths
        w, r, b = prices(prim, ins, res; Kg = K0, Lg = L0)
    #    println("Prices are\n wage = $w,\n rate $r, \n and benefits are $b")
        ugrid_work = utility_grid_work(prim, ins, res)
        ugrid_ret = utility_grid_ret(prim, ins, res)


        for t = (T-1):-1:1 ### do not include the last period since it is already known
        #    println("period t = $t")
            wage = w[t]
            rate = r[t]
            ben = b[t]

            backward_induct(prim, ins, res, t, ugrid_work, ugrid_ret)
        end

        # Find stationary distribution
        get_ss_distr(prim, ins, res, res.Γ[1, :, :, :])
        # Find aggregate variables according to the obtained decision rule and distribution
        K1, L1 = CalcAggregate(prim, ins, res)

        # Compare
        max_dif_K = maximum(abs.(K1 .- K0))
        max_dif_L = maximum(abs.(L1 .- L0))
        println("Max dif K is $max_dif_K and max dif L is $max_dif_L") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")


        dif = max_dif_K + max_dif_L
        println("FINDS DIFFERENCE $dif") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")
    #    println("K0 = $K0,\n K1 = $K1, \n L0 = $L0, \n L1 = $L1")

        ### We should not change the starting point!
        K1[1] = K0[1]
        L1[1] = L0[1]

        if dif > tol
            K0 = λ*K0 + (1-λ)*K1
            L0 = λ*L0 + (1-λ)*L1
        else
            K = K1
            L = L1
            println("DONE!\n")
            break
        end

        n_iter += 1
    end

    ##############

    return ins, res, K, L     # what we need in order to make graphs
end
