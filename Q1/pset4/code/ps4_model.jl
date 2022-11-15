### more functions here ###


# Finds the optimal labor supply
function getLabor(a_index, z_index, ap_index, prim::primitives, res::TransitionResults, t::Int64, age::Int64)
    @unpack e, γ, R, δ, θ, a_grid  = prim
    @unpack w, r, b, θ_path = res

    l = (γ*(1-θ_path[t])*e[age, z_index]*w[t] - (1-γ)*((1+r[t])*a_grid[a_index] - a_grid[ap_index])) / ((1-θ_path[t])*w[t]*e[age, z_index])

    if l > 1.0
        l = 1.0
    elseif l < 0.0
        l = 0.0
    end

    return l
end

#utility grid for retirees
function utility_grid_ret(prim::primitives, res::TransitionResults)
    @unpack e, γ, R, δ, σ, a_grid, T, na  = prim
    @unpack r, b = res

    ugrid_ret = zeros(na, na, T) #one for each a and a' and t

    for t = 1:T, a_index=1:na, ap_index=1:na
        a, ap = a_grid[a_index], a_grid[ap_index]
        c = (1+r[t])*a + b[t] - ap
        if c>0
            ugrid_ret[a_index, ap_index, t] = c^((1-σ)*γ)/(1-σ)
        else
            ugrid_ret[a_index, ap_index, t] = -Inf
        end
    end
    ugrid_ret
end

#utility grid for workers
function utility_grid_work(prim::primitives,  res::TransitionResults)
    @unpack N, R, σ, β, π, γ, e, z, a_grid, na, t_noss, θ, T  = prim
    @unpack w, r, b, θ_path = res

    ugrid_work = zeros(na*2, na, R, T) #current a and z x choice of ap x age x time
    # 4 dimension u grid, to keep track of age until retirement (R) and time

    #loop over time, current a and z, age, and  a'
    for t = 1:T, a_index = 1:na, age = 1:(R-1), ap_index = 1:na, z_index = 1:2
        a = a_grid[a_index]
        ap = a_grid[ap_index]
        labor = getLabor(a_index, z_index, ap_index, prim, res, t, age) #use optimal labor function
        c = w[t]*(1-θ_path[t])*e[age, z_index]*labor+(1+r[t])*a - ap

        if c>0
            inside=(c^γ)*(1-labor)^(1-γ)
            ugrid_work[a_index+na*(z_index - 1), ap_index, age, t] = (inside^(1-σ))/(1-σ)
        else
            ugrid_work[a_index + na*(z_index - 1), ap_index, age, t] = -Inf
        end
    end
    ugrid_work #deliverable
end

# Bellman function for retirees - let make age an argument of the function - maybe it will make the code faster
function bellman_retiree(prim::primitives, res::TransitionResults, t::Int64, age::Int64, ugrid_ret::Array{Float64,3})
    @unpack N, R, σ, β, π, γ, e, z, a_grid, na, t_noss, θ  = prim
    @unpack w, r, b, pol_func, val_func = res


    value_temp = zeros(na, 2)      # value function guess


        for z_index = 1:2, a_index = 1:na
            max_val = -Inf
            if age == N
                max_val = ugrid_ret[a_index, 1, t]
            else
                for ap_index in 1:na
                    val = ugrid_ret[a_index, ap_index, t] + β*val_func[t+1, ap_index, z_index, age+1]
                    if val > max_val
                        max_val = val
                        pol_func[t, a_index, z_index, age] = a_grid[ap_index]
                    end
                end
            end
            value_temp[a_index, z_index] = max_val
        end
    return value_temp
end

# Bellman function for workers
function bellman_worker(prim::primitives, res::TransitionResults, t::Int64, age::Int64, ugrid_work::Array{Float64,4})
    @unpack N, R, σ, β, π, γ, e, z, a_grid, na, T, t_noss, θ  = prim
    @unpack w, r, b, pol_func, labor, val_func = res


    value_temp = zeros(na, 2)      # value function guess

    for z_index = 1:2, a_index = 1:na
        max_val = -Inf

        for ap_index in 1:na
            val = ugrid_work[a_index + na*(z_index-1), ap_index, age, t]+β*sum(π[z_index,:].*val_func[t+1, ap_index, :, age+1])
            if val > max_val
                max_val = val
                pol_func[t, a_index, z_index, age] = a_grid[ap_index]
                labor[t, a_index, z_index, age] = getLabor(a_index, z_index, ap_index, prim, res, t, age)
            end
        end
        value_temp[a_index, z_index] = max_val
        end

    return value_temp
end

#backward induction protocol
function backward_induct(prim::primitives, res::TransitionResults, t::Int64, ugrid_work::Array{Float64,4}, ugrid_ret::Array{Float64,3})
    @unpack R, N, T = prim

    #loop over all ages
    for j = 1:N
        age = N-j+1 #now going backwards = 1 when j = J
        if age>= R
            res.val_func[t,:,:,age] = bellman_retiree(prim, res, t, age, ugrid_ret)
        else
            res.val_func[t,:,:,age] = bellman_worker(prim, res, t, age, ugrid_work)
        end
    end
end

# Create population distribution by age
function AgeDistribution(prim::primitives)
    @unpack N, n = prim

    μ = ones(N)
    μ[1] = 1

    for i = 2:N
        μ[i] = μ[i-1]/(1+n)
    end
    μ = μ ./sum(μ)
    return μ
end

# Generate stationary distribution
function get_ss_distr(prim::primitives, res::TransitionResults, initial_distr::Array{Any, 3})
    @unpack π, π0, N, n, na, a_grid, T = prim
    @unpack pol_func = res

    distr = zeros(T, na, 2, N)
    distr[1, :, :, :] = initial_distr

    for t = 1:(T-1)
        distr[t, 1, :, 1] = prim.μ_1 .* π0  # we know what they should start with at age = 1
        for age = 1:(N-1)
            for z_index = 1:2, a_index = 1:na
                    ap = pol_func[t, a_index, z_index, age]
                    for ap_index = 1:na, zp_index = 1:2 #loop over tomorrow's states
                        if a_grid[ap_index] == ap
                            distr[t+1, ap_index, zp_index, age+1] += distr[t, a_index, z_index, age]*
                                π[z_index, zp_index]*(1/(1+n))
                        end
                    end
            end
        end
    end
    res.Γ = distr
end

function CalcAggregate(prim::primitives, res::TransitionResults)
     @unpack R, N, e,  T, na, a_grid  = prim
     @unpack Γ, labor, K, L = res


     K1 = zeros(T)
     L1 = zeros(T)

     K1[1] = K[1]
     L1[1] = L[1]

     # Calculate aggregate capital demand (K1)
     # Calculate aggregate labor demand (L1)
    for t = 2:T
        for age = 1:N
            for z_index = 1:2, a_index = 1:na
                K1[t] += Γ[t, a_index, z_index, age]*a_grid[a_index]
                if age < R
                    L1[t] += Γ[t, a_index, z_index, age] * e[age, z_index] * labor[t, a_index, z_index, age]
                end
            end
        end
    end

    return K1, L1
end

function update_prices(prim::primitives, res::TransitionResults; K0::Array{Float64}, L0::Array{Float64}, μ::Array{Float64})
    @unpack α, δ, R, N, n, θ, T, t_noss = prim
    #@unpack K, L  = res

    # Path of wages w
    for t = 1:T
        res.w[t] = (1-α) * (K0[t]^α) * (L0[t]^(-α))

        # Path of interest rates r
        res.r[t] = α * (K0[t]^(α-1))*(L0[t]^(1-α)) - δ
    end

    # Path of pensions
    for t = 1:T
        if t < t_noss
            @unpack w = res
            res.θ_path[t] = θ
            res.b[t] = (θ*w[t]*L0[t])/ sum(μ[R:N])
        else
            res.θ_path[t] = 0
            res.b[t] = 0
        end
    end
    return res.w, res.r, res.b, res.θ_path
end

# Calculate EV
function CalcEV(prim::primitives, res::TransitionResults, res_old::StationaryResults, μ; al = 0.0, au = 20.0, na = 200, θ = 0.11, z = [3.0, 0.5], γ = 0.42, T = 50, t_noss = 2)
    @unpack N = prim

    # Initialize
    EV = zeros(prim.na, 2, N)
    EV_age = zeros(N)
    voters = zeros(N)

    for age = 1:N
        for z_index = 1:2, a_index = 1:na
            EV_numerator = res.val_func[1, a_index, z_index, age]
            EV_denominator = res_old.val_func[a_index, z_index, age]
            EV[a_index, z_index, age] = (EV_numerator/EV_denominator)^(1/prim.γ*(1-prim.σ)) - 1
        
            EV_age[age,1] += EV[a_index, z_index, age] * (res_old.Γ[a_index, z_index, age]/μ[age])

            voters[age,1] += (EV[a_index, z_index, age] >= 0) * (res_old.Γ[a_index, z_index, age]/μ[age])
        end
    end

   EV_age, voters
end

# Solve model
function solveModel(res_old::StationaryResults, res_new::StationaryResults; al = 0.0, au = 20.0, na = 200, θ = 0.11, z = [3.0, 0.5], γ = 0.42, T = 30, t_noss = 2, tol = 0.005, iter = 1000, λ = 0.5)
    prim, res = InitializeTrans(al, au, na, θ, z, γ, T, t_noss)

    @unpack N, R, σ, β, π, γ, e, z, a_grid, na, T, t_noss, θ = prim
    @unpack labor, pol_func, val_func, K, L, w, r, b, θ_path, Γ = res

    # Set end points based on stationary results
    K[1] = res_old.K
    K[T] = res_new.K
    L[1] = res_old.L
    L[T] = res_new.L

    # Set initial transition path for aggregate K and L
    K_step = (K[T] - K[1])/(T-1)
    L_step = (L[T] - L[1])/(T-1)
    for i = 2:T-1
        K[i] = K[i-1] + K_step
        L[i] = L[i-1] + L_step
    end

    # Population shares
    μ = AgeDistribution(prim)
#    println("mu = $μ")

    # Find prices
    w, r, b, θ_path = update_prices(prim, res; K0 = K, L0 = L, μ = μ)

    #fill in the terminal value function, stat distributions, policy functions
    res.val_func[T, :,:,:] = res_new.val_func
    res.pol_func[T, :,:,:] = res_new.pol_func
    res.labor[T, :,:,:] = res_new.labor
    res.Γ[T,:,:,:] = res_new.Γ

    diff = 100
    n_iter = 1

    # Solve the agent's problem
    while (diff > tol && n_iter < iter)
        println("Iter = $n_iter")

        #println("Prices are\n wage = $w,\n rate $r, \n and benefits are $b")
        ugrid_work = utility_grid_work(prim, res)
        ugrid_ret = utility_grid_ret(prim, res)
        #println("grids are filled in\n")

        #
        for t = (T-1):-1:1
            backward_induct(prim, res, t, ugrid_work, ugrid_ret)
        end

        # Find stationary distribution
        get_ss_distr(prim, res, res_old.Γ)
        # Find aggregate variables according to the obtained decision rule and distribution
        K1, L1 = CalcAggregate(prim, res)
        println("K0 = $K and\n L0 = $L")
        println("K1 = $K1 and\n L1 = $L1")

        # Compare
        max_diff_K = maximum(abs.(K1 .- K))
        max_diff_L = maximum(abs.(L1 .- L))
        println("Max diff K is $max_diff_K and max diff L is $max_diff_L") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")

        diff = max_diff_K + max_diff_L
        println("FINDS DIFFERENCE $diff") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")
    #    println("K0 = $K0,\n K1 = $K1, \n L0 = $L0, \n L1 = $L1")

        ### We should not change the starting point!
        K1[1] = K[1]
        L1[1] = L[1]

        if diff > tol
            K = (1-λ)*K1 + λ*K
            L = (1-λ)*L1 + λ*L
        end

        res.w, res.r, res.b, res.θ_path = update_prices(prim, res; K0 = K, L0 =  L, μ = μ)
        n_iter += 1

        if diff < tol
            res.K = K1
            res.L = L1
            println("DONE!\n")
            break
        end
    end

    EV_age, voters = CalcEV(prim, res, res_old, μ)

    ##############

    return prim, res, K, L, EV_age, voters     # what we need in order to make graphs
end