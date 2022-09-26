# Anna: Not sure whether I need to include some files here. Probably, not.

@unpack pol_func, val_func = res
@unpack a_grid = par

# exercise 1.
# savings of the worker of 20 years old with high ans low productivity shocks
# Anna: I don't like the policy function
StatsPlots.plot(a_grid, pol_func[:, 1, 20], label = "High productivity", xlabel = "Current assets, a", ylabel = "Future assets, a", title = "Savings of the worker at 20 years old")
StatsPlots.plot!(a_grid, pol_func[:, 2, 20], label = "Low productivity", legend =:bottomright)
savefig(figpath*"savings_at_20.png")

# value function of a retired of 50 years old
StatsPlots.plot(a_grid, val_func[:, 1, 50], xlabel = "Current assets, a", ylabel = "Value", title = "Value function of the retired at the age 50")
savefig(figpath*"value_at_50.png")