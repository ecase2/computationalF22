### Lorenz curve
function Lorenz(prim::Primitives, res::Results)
    @unpack pol_func, distr = res
    @unpack s_grid, na, ns = prim

    prob = zeros(na)
    cdf = zeros(na)
    wealth = zeros(na, ns)
    wealth[:, 1] = pol_func[:, 1] .+ 1
    wealth[:, 2] = pol_func[:, 2] .+ 0.5

    total_wealth = zeros(na, ns)
    total_wealth[:, 1] = distr[:, 1].*wealth[:, 1]
    total_wealth[:, 2] = distr[:, 2].*wealth[:, 2]
    aggregate_wealth = sum(total_wealth)
    total_wealth = total_wealth/aggregate_wealth

    N = na*ns
    wealth_stack = zeros(N)
    wealth_stack = [total_wealth[:, 1]; total_wealth[:, 2]]

    people_sh = [distr[:, 1]; distr[:, 2]]

    df = DataFrame([wealth_stack, people_sh], :auto)
    rename!(df, [:x1, :x2] .=>  [:wealth, :people_sh])
    df = sort!(df) ### Important
    sum_people = zeros(N)
    sum_wealth = zeros(N)

    for i in 1:N
        sum_people[i] = sum(df.people_sh[1:i])
        sum_wealth[i] =  sum(df.wealth[1:i])
    end
    sum_people[N] = 1 ### Without it get 0.9999999999999887
    sum_wealth[N] = 1 ### Without it get 0.950845200818719 since negative wealth is also positive
    return sum_wealth, sum_people, df, wealth_stack
end

function Gini(prim::Primitives, res::Results; wealth = sum_wealth, share = sum_people)
    @unpack na, ns = prim

    ###
    area0 = 0.5
    area1 = 0
    x = collect(0:1/(ns*na-1):1)
    y = x
    diff = zeros(ns*na)
    for i in 1:length(share)
        if wealth[i] > 0
            diff[i] = y[i] - wealth[i]
        end
        if wealth[i] <= 0
            diff[i] = wealth[i] + y[i]
        end
    end
    for i in 1:(length(share)-1)
        area1 = area1 + (share[i+1] - share[i])*(diff[i] + diff[i+1])
    end

    return area1, area1/area0

end

###### To plot the Lorenz curve and compute Gini
sum_wealth, sum_people, df, wealth_stack = Lorenz(prim, res)


x = collect(0:1/(prim.ns*prim.na-1):1)
Plots.plot(x, x, title="Lorenz curve", label = "45 degree line")
Plots.plot!(sum_people, sum_wealth, title="Lorenz curve", label = "Lorenz curve", legend =:bottomright)

area, G = Gini(prim, res; wealth = sum_wealth, share = sum_people)
