#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 5
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

# Initialize
function Initialize()
    par = Params()      # initialize parameters
    grid = Grids()      # initialize grids
    shocks = Shocks()   # initialize shocks draw

    # Set guesses
    pf_k = zeros(grid.n_k, grid.n_eps, grid.n_K, grid.n_z)
    pf_v = zeros(grid.n_k, grid.n_eps, grid.n_K, grid.n_z)

    a0 = 0.095
    a1 = 0.999
    b0 = 0.085
    b1 = 0.999

    R2 = 0.0

    res = Results(pf_k, pf_v, a0, a1, b0, b1, R2)   # initialize results structure

    return par, grid, shocks, res
end

# Calculate prices (wage, interest rate)
function calc_prices(par::Params, z::Float64, K::Float64, L::Float64)
    @unpack cALPHA = par

    # Calculate wage 
    w = (1-cALPHA) * z * K^cALPHA * L^(-cALPHA)

    # Calculate interest rate 
    r = cALPHA * z * K^(cALPHA-1) * L^(1-cALPHA)

    w, r
end

# Value function iteration
function V_iterate(par::Params, grid::Grids, shocks::Shocks, res::Results; iter = 1000, tol = 0.005)
    n = 0 # counter
    diff = 10.0

    while diff > tol 
        pf_k_next, pf_v_next = Bellman(par, grid, shocks, res)
        diff = maximum(abs.(pf_v_next .- res.pf_v))

        res.pf_v = pf_v_next # update value function
        res.pf_k = pf_k_next # update policy function

        n+=1
    end
end

# Simulate aggregate capital path
function simulate_KPath(par::Params, grid::Grids, res::Results, E::Array{Float64,2}, Z::Array{Float64,1}; K_ss = 11.55, drop_rows = 1000)
    @unpack N, T = par
    @unpack K_grid, z_grid, k_grid, eps_grid = grid
    @unpack pf_k = res

    # Initialize starting values of capital path and the capital matrix using steady state value
    KMatrix = zeros(N,T)
    KMatrix[:,1] .= K_ss

    KPath = zeros(T)
    KPath[1] = K_ss

    for t in 2:T
        z = Z[t]

        for n in 1:N
            eps = E[n,t]

            k_z_eps = pf_k[:, Integer(eps), :, Integer(z)]
            k_interp = interpolate(k_z_eps, BSpline(Linear()))

            k_index = get_index(KMatrix[n,t-1], k_grid)
            K_index = get_index(KPath[t-1], K_grid)
            
            KMatrix[n,t] = k_interp[k_index, K_index]
        end

        KPath[t] = sum(KMatrix[:, t])/N  
    end

    KPath = hcat(KPath, Z)
    KPath = KPath[drop_rows+1:T, :]

    return KPath
end

# Create matrix of aggregate K yesterday and aggregate K today
function reshape_K(K::Array{Float64,2}; K_ss = 11.55)
    A = zeros(length(K[:,1]))

    A[1] = K_ss
    for t = 2:length(K[:,1])
        A[t] = K[t-1]
    end

    B = hcat(A, K)

    C = B[sortperm(B[:,3]),:]

    return C
end

# Estimate regression
function estimate_reg(res::Results, K::Array{Float64,2})

    K = reshape_K(K)

    # Identify regression for good shocks
    y_a = log.(K[ K[:,3] .== 2.0, 1])
    x_a = [ones(length(y_a)) log.(K[ K[:,3] .== 2.0, 2])]

    # Identify regression for bad shocks
    y_b = log.(K[ K[:,3] .== 1.0, 1])
    x_b = [ones(length(y_b)) log.(K[ K[:,3] .== 1.0, 2])]

    # Find estimates
    beta_a = inv(x_a' * x_a) * x_a' * y_a
    beta_b = inv(x_b' * x_b) * x_b' * y_b

    ssr_a = sum((y_a - x_a * beta_a).^2)
    sst_a = sum((y_a .- mean(y_a)).^2)
    ssr_b = sum((y_b - x_b * beta_b).^2)
    sst_b = sum((y_b .- mean(y_b)).^2)

    res.R2 = 1 - (ssr_a + ssr_b)/(sst_a + sst_b)

    return beta_a[1], beta_a[2], beta_b[1], beta_b[2]
end

# Solve model
function solve_model(; λ = 0.5, iter = 1000, tol = 0.005)    
    # Initialize
    par, grid, shocks, res = Initialize()
    @unpack pf_k, pf_v = res
    
    # Draw shocks
    E, Z = draw_shocks(shocks, par.N, par.T)

    # Set initial guess
    a0_old = 0.095
    a1_old = 0.999
    b0_old = 0.085
    b1_old = 0.999

    diff = 10.0
    n = 1

    # Iterate
    while (diff > tol && n < iter)
        println("BEGINNING ITERATION $n")
        n = n+1

        # Solve
        V_iterate(par, grid, shocks, res)

        # Simulate capital path
        K = simulate_KPath(par, grid, res, E, Z)

        # Estimate regression
        a0_new, a1_new, b0_new, b1_new = estimate_reg(res, K)

        # Calculate difference
        diff = abs(a0_new - a0_old) + abs(a1_new - a1_old) + abs(b0_new - b0_old) + abs(b1_new - b1_old)
        println("Difference: $diff")
        println("\ta0_old = $a0_old \t a0_new = $a0_new")
        println("\ta1_old = $a1_old \t a1_new = $a1_new\n")
        println("\tb0_old = $b0_old \t b0_new = $b0_new")
        println("\tb1_old = $b1_old \t b1_new = $b1_new\n")

        # Adjust guess
        if diff > tol
            a0_old = λ*a0_new + (1-λ)*a0_old
            a1_old = λ*a1_new + (1-λ)*a1_old
            b0_old = λ*b0_new + (1-λ)*b0_old
            b1_old = λ*b1_new + (1-λ)*b1_old
        else
            println("DONE!\n")
            break
        end
    end

    output = [a0_old; a1_old; b0_old; b1_old; res.R2]

    return output
end

# Write output table
function write_table(output)
    set_default(fmt = "%.3f",           # latexify 
    convert_unicode = false, 
    latex = false)  

    # define the latex table head and foot 
    tabhead = "
    \\begin{table}\\caption{Final parameter estimates and \$R^2\$}\n\\centering
    \\begin{tabular}{ccccc}
    \\toprule
    \t \$a_0\$ & \$a_1\$ & \$b_0\$ & \$b_1\$ & \$R^2\$ \\\\
    \\hline"

    tabfoot = "
    \\bottomrule
    \\end{tabular}
    \\end{table}"

    # use latexify to create the table part (minus the column names, those are too complicated for latexify) 
    tabvals   = output
    temptable = latexify(tabvals, env = :tabular)

    # get rid of latexify's table head and foot, so that we can replace it with our own. 
    temptable = temptable[25:end-14]
    # put everything together and write the table 
    write(joinpath(root, "ps5_results.tex"), tabhead*temptable*tabfoot)
end