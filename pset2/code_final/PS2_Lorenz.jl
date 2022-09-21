### Lorenz curve
function Lorenz(prim::Primitives, res::Results)
    @unpack pol_func, distr = res
    @unpack s_grid, a_grid,  na, ns = prim

    N = na*ns

    wealth = zeros(na, ns)
    wealth[:, 1] = a_grid .+ 1
    wealth[:, 2] = a_grid .+ 0.5

    wealth = sort([wealth[:, 1]; wealth[:, 2]])

    people_sh = zeros(N)

    for w_index =  1:N, a_index = 1:na
        for s_index in 1:2
            w = wealth[w_index]
            if a_grid[a_index] + s_grid[s_index] == w
                people_sh[w_index] = people_sh[w_index] + distr[a_index, s_index]
            end
        end
    end

    total_wealth = zeros(N)
    total_wealth = people_sh.*wealth

    agg_wealth = sum(total_wealth)

    sum_people = zeros(N)
    sum_wealth = zeros(N)

    for i in 1:N
        sum_people[i] = sum(people_sh[1:i])
        sum_wealth[i] =  sum(total_wealth[1:i])/agg_wealth
    end
    return sum_wealth, sum_people
end


function Gini(prim::Primitives, res::Results; wealth = sum_wealth, share = sum_people)
    @unpack na, ns = prim

    ### Area under the 45 degree line
    area0 = 0.5
    ### Area under the Lorenz curve
    area1 = 0

    dif = zeros(ns*na)
    for i in 1:length(sum_wealth)
            dif[i] = sum_people[i] - wealth[i]
    end
    for i in 1:(length(sum_wealth)-1)
        ### compute area as a sum of rectangular areas
        area1 = area1 + (share[i+1] - share[i])*dif[i]

    end

    return area1/area0
end
