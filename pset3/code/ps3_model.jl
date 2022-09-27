# construct the grid for consumption and labor
# it should make the code faster
function fill_end_grids(par::parameters, res::results, grid::grids)
    @unpack a_grid, na, z, nz, N, R, θ, γ = par
    @unpack e, w, r, b = res
    @unpack c, l = grid

    for a_index = 1:na, z_index = 1:nz, j = 1:N, ap_index = 1:na
        if j < R
            grid.l[a_index, z_index, j, ap_index] = (γ*(1-θ)*e[j, z_index]*w - (1-γ)*((1+r)*a_grid[a_index] - a_grid[ap_index])) / ((1-θ)*w*e[j, z_index])
            if grid.l[a_index, z_index, j, ap_index] < 0
                grid.l[a_index, z_index, j, ap_index] = 0
            end
            if grid.l[a_index, z_index, j, ap_index] > 1
                grid.l[a_index, z_index, j, ap_index] = 1
            end
            grid.c[a_index, z_index, j, ap_index] = w*(1-θ)*e[j, z_index]*l[a_index, z_index, j, ap_index] + (1+r)*a_grid[a_index] - a_grid[ap_index]
            # important: le may be negative according to the formula - we will need to proceed only if l is between 0 and 1.
        end
        if j >= R
            #grid.l[a_index, z_index, j, ap_index] = 0 # redundant
            grid.c[a_index, z_index, j, ap_index] = (1+r)*a_grid[a_index] + b - a_grid[ap_index]
        end
    end
end


function backward_iteration(par::parameters, res::results, grid::grids)
    @unpack a_grid, na, z, nz, N, R, θ, γ, σ, β, π = par
    @unpack e, w, r, b, val_func, pol_func = res
    @unpack c, l = grid

    # last period - save nothing, consume everythinig
    for a_index = 1:na, z_index = 1:nz
        pol_func[a_index, z_index, N] = 0
        cons                          = c[a_index, z_index, N, 1]
        val_func[a_index, z_index, N] = cons^((1-σ)*γ)/(1-σ)
    end
    t = N-1
    while t > 0
        println("age is ", t)
        if t >= R
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
        if t < R
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

function get_distr(res::results, par::parameters)
    @unpack pol_func = res 
    @unpack π, na, a_grid, N, n = par 
    F = zeros(par.na, par.nz, par.N)
    F[1,:, 1] = par.μ_1 .* [0.2037 0.7963]  # we know what they should start with at age = 1
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

## INSERT MARKET CLEARING FUNCTION HERE

function solve(;modeltype::String = "benchmark")
    par, res, grid = Initialize(modeltype)
    fill_end_grids(par, res, grid)
    backward_iteration(par, res, grid)
    get_distr(res, par)

    # USE MARKET CLEARING FUNCTION HERE 

    # NEED TO OUTPUT: K, L, w, r, b, W, cv 
end