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

    R2 = zeros(2)

    res = Results(pf_k, pf_v, a0, a1, b0, b1, R2)   # initialize results structure

    par, grid, shocks, res
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

# Simulate aggregate capital path
function simulate_KPath(par::Params, grid::Grids, res::Results; K_ss = 11.55)
    @unpack N, T = par
    @unpack K_grid, z_grid, k_grid, eps_grid = grid
    @unpack pf_k = res

    KMatrix = zeros(N,T+1)
    KMatrix[:,1] .= K_ss

    KPath = zeros(T+1)
    KPath[1,1] = K_ss

    for t in 2:T+1
        KBar = KPath[t-1,1]
        K_index = get_index(KBar, K_grid)

        z = Z[t,1]
        z_index = get_index(z, z_grid)
        
        for n in 1:N
            k = KMatrix[n,t-1]
            k_index = get_index(k, k_grid)

            eps = E[n,t]
            eps_index = get_index(eps, eps_grid)
            
            KMatrix[n,t] = pf_k[k_index, eps_index, K_index, z_index]
        end

        KPath[t,1] = sum(A, dims = 1)  
    end

    
        

    end

    



    V = zeros(N)

    for t in 1:T
        K = zeros(N)

        pf_k = zeros(grid.n_k, grid.n_eps, grid.n_K, grid.n_z)



        V = hcat(V,K)
    end

    V = V[N,2:T]

    V
end

# Estimate regression
function estimate_reg()
    # Sort based on z

    # Estimate regression

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
        pf_k, pf_v = Bellman(par, grid, shocks, res)

        # Simulate capital path


        # Estimate regression
        a0_new, a1_new, b0_new, b1_new = estimate_reg()

        # Calculate difference
        diff = abs(a0_new - a0_old) + abs(a1_new - a1_old) + abs(b0_new - b0_old) + abs(b1_new - b1_old)

        # Adjust guess
        if diff > tol
            a0_old = λ*a0_new + (1-λ)*a0_old
            a1_old = λ*a1_new + (1-λ)*a1_old
            b0_old = λ*b0_new + (1-λ)*b0_old
            b1_old = λ*b1_new + (1-λ)*b1_old
        else
            #[Insert something here]

            println("DONE!\n")
            break
        end
    end

    return 
end