#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 3
    AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

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

function getConsumption(a_index, z_index, age::Int64, ap_index, par::parameters, res::results)
    @unpack w, θ, e, r, b, γ = res
    @unpack a_grid, R, δ = par
    if age < R
        l = getLabor(a_index, z_index, age, ap_index, par, res)
        c = w*(1-θ)*e[age, z_index]*l + (1+r-δ)*a_grid[a_index] - a_grid[ap_index] # γ * ( (1 - θ)*e[age, z_index]*w + (1+r)*a_grid[a_index] - a_grid[ap_index])# 
        return c, l
    else
        c = (1+r-δ)*a_grid[a_index] + b - a_grid[ap_index]
        return c
    end

end

# Do backwards induction to get value function and policy function
function backward_iteration(par::parameters, res::results)
    @unpack a_grid, na, nz, N, R, σ, β, π = par
    @unpack val_func, pol_func, θ, γ, e, z = res # w, r, b,

    # last period - save nothing, consume everythinig
    for a_index = 1:na, z_index = 1:nz
        pol_func[a_index, z_index, N] = 0.0
        cons                          = getConsumption(a_index, z_index, N, 1, par, res)
        if cons > 0.0
            val_func[a_index, z_index, N] = cons^((1-σ)*γ)/(1-σ)
        end
    end

    # For retirees
    for age = N-1:-1:R
        #println("age is ", age)
        for a_index = 1:na, z_index = 1:nz
            max_val = -Inf
            for ap_index = 1:na
                cons = getConsumption(a_index, z_index, age, ap_index, par, res)
                if cons >= 0
                    util = cons^((1-σ)*γ)/(1-σ) + β*val_func[ap_index, z_index, age + 1]
                    if util > max_val
                        max_val = util
                        pol_func[a_index, z_index, age] = a_grid[ap_index]
                        val_func[a_index, z_index, age] = max_val
                    end 
                end 
            end 
        end 
    end 

    # For workers
    for age = R-1:-1:1
        #println("age is ", age)
        for a_index = 1:na, z_index = 1:nz
            max_val = -Inf
            for ap_index = 1:na
                cons, lab = getConsumption(a_index, z_index, age, ap_index, par, res)
                if (cons >= 0) && (0<= lab <= 1)
                    util = (cons^γ*(1-lab)^(1-γ))^(1-σ)/(1-σ) + β*(transpose(π[z_index, :])*val_func[ap_index, :, age + 1])
                    if util > max_val
                        max_val = util
                        pol_func[a_index, z_index, age] = a_grid[ap_index]
                        res.l[a_index, z_index, age] = lab
                        val_func[a_index, z_index, age] = max_val
                    end
                end
            end
        end
    end
end




# utility_retiree: this function encodes the utility function of the retiree
function utility_retiree(c::Float64, σ::Float64, γ::Float64)
    if c > 0
        u_r = (c^((1-σ)*γ))/(1-σ) # CRRA utility for retiree
    else
        u_r = -Inf
    end
    u_r # return calculation
end

function utility_worker(c::Float64, l::Float64, σ::Float64, γ::Float64)
    if c > 0
        u_w = (((c^γ)*(1-l)^(1-γ))^(1-σ))/(1-σ) # CRRA utility for worker
    else
        u_w = -Inf
    end
    u_w # return calculation
end

# bellman_retiree: this function encodes the Bellman function for the retired
function bellman_retiree(prim::Primitives, res::Results, age::Int64)
    @unpack a_grid, na, nz, N, R, σ, β, π = par
    @unpack val_func, pol_func, θ, γ, e, z = res

    val_index = retiree_val_index(age_retire, nz, age)                  # mapping to index in val func

    if age == N # if age == N, initialize last year of life value function
        for (a_index, a_today) in enumerate(a_grid)
            c = (1 + r) * a_today + b                                   # consumption in last year of life (a' = 0)
            res.val_func[a_index, val_index] = utility_retiree(c, σ, γ) # value function for retiree given utility (v_next = 0 for last period of life)
        end
    else # if not at end of life, compute value funaction for normal retiree
        choice_lower = 1                                                # for exploiting monotonicity of policy function

        for (a_index, a_today) in enumerate(a_grid)                     # loop through asset levels today
            max_val = -Inf                                              # initialize max val

            @sync @distributed for ap_index in choice_lower:na          # loop through asset levels tomorrow
                a_tomorrow = a_grid[ap_index]                           # get a tomorrow
                v_next_val = res.val_func[ap_index, val_index+1]        # get next period val func given a'
                c = (1 + r) * a_today + b - a_tomorrow                  # consumption for retiree

                if c > 0                                                    # check for positivity of c
                    v_today = utility_retiree(c, σ, γ) + β * v_next_val     # value function for retiree

                    if max_val < v_today  # check if we have bigger value for this a_tomorrow
                        max_val = v_today                                   # update max value
                        res.pol_func[a_index, val_index] = a_tomorrow       # update asset policy
                        choice_lower = ap_index
                    end
                end
            end # end of a_tomorrow loop for standard retiree
            res.val_func[a_index, val_index] = max_val # update val function for a_today and age
        end # end of a_today loop for standard retiree
    end # end of if statement checking for whether at end of life or not
end

function bellman_worker(prim::Primitives, res::Results, age::Int64)
    @unpack a_grid, na, nz, N, R, σ, β, π = par
    @unpack val_func, pol_func, θ, γ, e, z = res

    for z_index = 1:nz

        choice_lower = 1

        for a_index = 1:na
            a = a_grid[a_index]
            maxval = -Inf

            for ap_index in choice_lower:na
                a_prime = a_grid[ap_index]                           
                l = labor_supply(γ, θ, e_today, w, r, a_today, a_tomorrow)       
                c = w * (1-θ) * e_today * l + (1 + r) * a - a_prime

                if c > 0 && l >= 0 && l <= 1                            # check for positivity of c and constraint on l: 0 <= l <= 1
                    if age == R-1
                        v_next_val = res.val_func[ap_index, age * nz + 1]       # get next period val func (just scalar) given a_tomorrow
                        v_today = utility_worker(c, l, σ, γ) + β * v_next_val
                    end
                end
            end

        end


    end



    for (z_index, z_today) in enumerate(z)                              # loop through productivity states
        e_today = e[age, z_index]                                       # get productivity for age and z_today
        z_prob = z_matrix[z_index, :]                                   # get transition probabilities given z_today
        val_index = worker_val_index(z_index, age, nz)                  # get index to val, pol, lab func

        choice_lower = 1                                                # for exploiting monotonicity of policy function

        for (a_index, a_today) in enumerate(a_grid)                     # loop through asset levels today
            max_val = -Inf                                              # initialize max val

            @sync @distributed for ap_index in choice_lower:na          # loop through asset levels tomorrow

                a_tomorrow = a_grid[ap_index]                           # get a tomorrow
                l = labor_supply(γ, θ, e_today, w, r, a_today, a_tomorrow)       # labor supply for worker
                c = w * (1-θ) * e_today * l + (1 + r) * a_today - a_tomorrow     # consumption for worker

                if c > 0 && l >= 0 && l <= 1                            # check for positivity of c and constraint on l: 0 <= l <= 1
                    if age == age_retire -1                                     # if age == 45 (retired next period), then no need for transition probs
                        v_next_val = res.val_func[ap_index, age * nz + 1]       # get next period val func (just scalar) given a_tomorrow
                        v_today = utility_worker(c, l, σ, γ) + β * v_next_val   # value function for worker

                    else                                                                  # else, need transition probs
                        v_next_val = res.val_func[ap_index, (age+1)*nz - 1:(age+1)*nz]    # get next period val func (vector including high and low prod) given a_tomorrow
                        v_today = utility_worker(c, l, σ, γ) + β * z_prob' * v_next_val   # value function for worker with transition probs
                    end

                    if max_val <= v_today  # check if we have bigger value for this a_tomorrow
                        max_val = v_today                                   # update max value
                        res.pol_func[a_index, val_index] = a_tomorrow       # update asset policy
                        res.lab_func[a_index, val_index] = l                # update labor policy
                        choice_lower = ap_index
                    end
                end
            end # end of a_tomorrow loop for worker
            res.val_func[a_index, val_index] = max_val # update val function for a_today and age
        end # end of a_today loop for worker
    end
end

# Backward iteration
function V_iterate(par::parameters, res::results)
    @unpack N, R = par

    for age in N:-1:1
        if age >= R                  # if between retirement age and end of life
            bellman_retiree(par, res, age)    # call Bellman for retiree
        else                                   # else, agent is still worker
            bellman_worker(par, res, age)     # call Bellman for worker
        end
    end
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
    @unpack α, δ, R, N = par
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
        res.b = (θ*w*L)/ sum(F[ :, :, R:N])
    end
end

# Determine market clearing by calculating aggregate capital demand and aggregate labor demand
function AggregateDemand(par::parameters, res::results)
    @unpack R, N, na, nz, a_grid = par
    @unpack F, l, e = res

    K1 = 0.0
    L1 = 0.0

    # Calculate aggregate capital demand (K1)
    K1 = sum(F .* a_grid)

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

# Calculate welfare objects - need to sum over result in order to calculate total welfare
function CalcWelfare(res::results)
    @unpack val_func, F = res

    # Calculate welfare
    welfare = val_func .* F
    welfare = sum(welfare)

    return welfare
end

# Calculate CV
function CalcCV(res::results)
    @unpack pol_func, F = res

    # Calculate wealth
    wealth = pol_func .* F

    # Calculate mean welfare
    wealth_mean = mean(wealth)

    # Calculate standard deviation
    wealth_sd = std(wealth)

    # Calculate CV
    cv = wealth_sd/wealth_mean

    return cv
end

function computecv(par::parameters, res::results)
    @unpack a_grid, na, nz, N, R = par
    @unpack θ, w, r, b, e, l, F = res

    mean_wealth = 0.0
    var_wealth  = 0.0
    cv_wealth   = 0.0

    for i=1:par.na, j=1:par.nz, age=1:par.N
        if age < par.R
            mean_wealth += (res.w * (1-res.θ) * e[age, j] * res.l[i,j,age] + (1+res.r) * par.a_grid[i]) * res.F[i,j,age]
        elseif age >= par.R
            mean_wealth += ((1 + res.r) * par.a_grid[i] + res.b) * res.F[i,j,age]
        end
    end


    for i=1:par.na, j=1:par.nz, age=1:par.N
        if age < par.R
            var_wealth += (res.w * (1 - res.θ) * e[age, j] * res.l[i,j,age] + (1+res.r) * par.a_grid[i] - mean_wealth)^2 * res.F[i,j,age]
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

    # Set (smart) initial guesses
    K0 = 4.7
    L0 = 0.4
    if sum(z) == 1.0
        K0 = 1.0
        L0 = 0.15
    elseif γ == 1.0
        K0 = 7.0
        L0 = 0.7
    end 

    diff = 10.0
    n = 1

    # Find K, L
    while (diff > tol && n < iter)
        println("BEGINNING ITERATION $n")
        n = n+1
        
        # Calculate w, r, b based on guess of K0, L0
        CalcPrices(par, res, K0, L0)
        w = res.w
        r = res.r
        b = res.b
        println("Wage = $w, Interest rate = $r, Benefit = $b")

        # Solve
        backward_iteration(par, res)
        get_distr(par, res)

        K1, L1 = AggregateDemand(par, res)
        diff = abs(K1 - K0) + abs(L1 - L0)
        println("FINDS DIFFERENCE $diff") # ... K DIFFERENCES $(K1 - K0), L DIFFERENCES $(L1 - L0)\n")
        println("\tK0 = $K0 \t K1 = $K1")
        println("\tL0 = $L0 \t L1 = $L1\n")


        # Adjust guess
        if diff > tol
            # if diff > 1.0
            #     if (K1 - K0) > 0 # if the demand is higher than supply, new supply guess has more weight
            #         λ_k = 0.3
            #     else
            #         λ_k = 0.7
            #     end
            #     if (L1 - L0) >0
            #         λ_l = 0.3
            #     else
            #         λ_l = 0.7
            #     end
            # else
                #λ = 0.5
            # end
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
    cv = CalcCV(res)
    cv2 = computecv(par, res)
    println("WELFARE IS $welfare ... CV IS $cv")

    # Produce results matrix
    output = [res.K; res.L; res.w; res.r; res.b; welfare; cv; cv2]

    if (θ == 0.11) && sum(z) >2.0 && (γ != 1.0 )
        return output, par, res
    else
        return output
    end
end
