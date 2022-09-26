# construct the grid for consumption and labor
# it should make the code faster
function fill_end_grids(par::parameters, res::results, grid::grids)
    @unpack a_grid, na, z, nz, N, R, θ, γ = par
    @unpack e, w, r, b = res
    @unpack c, l = grid

    for a_index = 1:na, z_index = 1:nz, j = 1:N, ap_index = 1:na
        if j < 46
            grid.l[a_index, z_index, j, ap_index] = (γ*(1-θ)*e[j, z_index]*w - (1-γ)*((1+r)*a_grid[a_index] - a_grid[ap_index]))/((1-θ)*w*e[j, z_index])
            if grid.l[a_index, z_index, j, ap_index] < 0
                grid.l[a_index, z_index, j, ap_index] = 0
            end
            if grid.l[a_index, z_index, j, ap_index] > 1
                grid.l[a_index, z_index, j, ap_index] = 1
            end
            grid.c[a_index, z_index, j, ap_index] = w*(1-θ)*e[j, z_index]*l[a_index, z_index, j, ap_index] + (1+r)*a_grid[a_index] - a_grid[ap_index]
            # important: le may be negative according to the formula - we will need to proceed only if l is between 0 and 1.
        end
        if j >= 46
            #grid.l[a_index, z_index, j, ap_index] = 0 # redundant
            grid.c[a_index, z_index, j, ap_index] = (1+r)*a_grid[a_index] + b - a_grid[ap_index]
        end
    end
end


function backward_iteration(par::parameters, res::results, grid::grids)
    @unpack a_grid, na, z, nz, N, R, θ, γ, σ, β, π = par
    @unpack e, w, r, b, val_func, pol_func, μ = res
    @unpack c, l = grid

    # last period - save nothing, consume everythinig
    for a_index = 1:na, z_index = 1:nz
        pol_func[a_index, z_index, N] = 0
        cons = c[a_index, z_index, N, 1]
        val_func[a_index, z_index, N] = cons^((1-σ)*γ)/(1-σ)
    end
    t = N-1
    while t > 0
        println("age is ", t)
        if t > 45
            for a_index = 1:na, z_index = 1:nz # not very efficient. Make to do nz loops which are identical
                max_val = -Inf
                for ap_index = 1:na
                    cons = c[a_index, z_index, t, ap_index]
                    if cons >= 0
                        util = cons^((1-σ)*γ)/(1-σ) + β*res.val_func[ap_index, z_index, t + 1]
                        if util > max_val
                            max_val = util
                            res.pol_func[a_index, z_index, t] = a_grid[ap_index]
                            res.val_func[a_index, z_index, t] = util
                        end
                    end
                end
            end
        end
        if t <= 45
            for a_index = 1:na, z_index = 1:nz
                max_val = -Inf
                for ap_index = 1:na
                    lab = l[a_index, z_index, t, ap_index]
                    cons = c[a_index, z_index, t, ap_index]
                    if cons >= 0
                        util = (cons^γ*(1-lab)^(1-γ))^(1-σ)/(1-σ) + β*(transpose(π[z_index, :])*res.val_func[ap_index, :, t + 1])
                        if util > max_val
                            max_val = util
                            res.pol_func[a_index, z_index, t] = a_grid[ap_index]
                            res.val_func[a_index, z_index, t] = util
                        end
                    end
                end
            end

        end
        t -=  1
    end

end

function rel_sizes(par::parameters, res::results, grid::grids)
    @unpack n, μ_1, N = par

    relsize = zeros(N)
    relsize[1] = μ_1
    for t = 2:N
        relsize[t] = relsize[t-1]/(1+n)
    end
    return relsize
end
# DO WE NEED RELATIVE SIZES TO COMPUTE STUFF? I think no, since we have found μ_1 using math

function transition_function(par::parameters, res::results, grid::grids)
    @unpack a_grid, na, z, nz, N, R, π = par
    @unpack val_func, pol_func = res

    ### Indicator for policy function
    χ = zeros(na, na, nz, N) # 1 - future assets, 2 - current assets, 3 - current productivity, 4 - current age
    for ap_index = 1:na, a_index = 1:na, z_index = 1:nz, t = 1:N
        if pol_func[a_index, z_index, t] == a_grid[ap_index]
            χ[ap_index, a_index, z_index, t] = 1
        end
    end
    # transition function is a product of the indicator function for policy rule and the transition matrix
    P = zeros(na, nz, na, nz, N) # 1 - future assets, 2 - future productivity, 3 - current assets, 4 - current productivity, 5 - current age
    for ap_index = 1:na, zp_index = 1:nz, a_index = 1:na, z_index = 1:nz, t = 1:N
        P[ap_index, zp_index, a_index, z_index, t] = χ[ap_index, a_index, z_index, t]*π[z_index, zp_index]
    end
    return P
end


function distribution(par::parameters, res::results, grid::grids; tol = 0.001, iter = 1000)
    @unpack a_grid, na, z, nz, N, R, θ, γ, σ, β, π, μ_1, pzᴴ, pzᴸ, n = par
    @unpack e, w, r, b, val_func, pol_func, μ = res
    @unpack c, l = grid

    relsize = rel_sizes(par, res, grid)
    P = transition_function(par, res, grid)

    # distribution of agents
    F = zeros(na, nz, N)

    # distribution of newborns - they all are born with 0 assets
    F[1, 1, 1] = pzᴴ*μ_1
    F[1, 2, 1] = pzᴸ*μ_1


    # using transition_function P
    for zp_index in 1:nz, ap_index = 1:na, t = 2:N
            F[ap_index, zp_index, t] = (transpose(P[ap_index, zp_index, :, 1, t-1])*(F[:, 1, t-1]) +
            transpose(P[ap_index, zp_index, :, 2, t-1])*(F[:, 2, t-1]))
    end

    # without transition_function P
    #    for a_index = 1:na, z_index = 1:nz, t = 2:N
    #        for ap_index = 1:na, zp_index = 1:nz
    #            if pol_func[a_index, z_index, t] == a_grid[ap_index]
    #                F[ap_index, zp_index, t] = F[ap_index, zp_index, t] + F[a_index, z_index, t-1]*π[z_index, zp_index]
    #            end
    #        end
    #    end

    F
end

# using the idea as in the pset2 - now sure why there we have used iteration until convergence,
# is it because  the model was inifinite horizon model?
# anyways, the code works wrong since it converges to distribution = zeros(na, nz, N)
function distribution_v2(par::parameters, res::results, grid::grids; tol = 0.001, iter = 1000)
    @unpack a_grid, na, z, nz, N, R, θ, γ, σ, β, π, μ_1, pzᴴ, pzᴸ, n = par
    @unpack e, w, r, b, val_func, pol_func, μ = res
    @unpack c, l = grid

    P = transition_function(par, res, grid)

    # distribution of agents
    F = zeros(na, nz, N)

    # distribution of newborns - they all are born with 0 assets
    F[1, 1, 1] = pzᴴ*μ_1
    F[1, 2, 1] = pzᴸ*μ_1

    F[2:na, :, 2:N] = (ones(na-1, nz, N-1) .- (F[1, 1, 1]+F[1, 2, 1]))/((na-1)*nz*(N-1))
    diff = 100
    n_it = 1
    ### create the indicator variable
    while (diff > tol && n_it < iter)
        println("iteration is ", n_it, " diff is ", diff)
        distr_1 = zeros(na, nz, N)
        for zp_index in 1:nz, ap_index = 1:na, t = 2:N
                distr_1[ap_index, zp_index, t] = (transpose(P[ap_index, zp_index, :, 1, t-1])*(F[:, 1, t-1]) +
                transpose(P[ap_index, zp_index, :, 2, t-1])*(F[:, 2, t-1]))
        end
        println("max(distr_1) = ", maximum(distr_1))
        diff = maximum(abs.(distr_1 .- F))
        n_it = n_it + 1
    #    println("Iteration ", n-1, " Diff = ", diff);
        F = distr_1
    end
    F
end
