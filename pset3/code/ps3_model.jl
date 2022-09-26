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
