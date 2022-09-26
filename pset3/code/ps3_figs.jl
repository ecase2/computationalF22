# set some defaults for the plots:
default(titlefont = (20, "times"), linewidth = 2)

# Anna: Not sure whether I need to include some files here. Probably, not.

@unpack pol_func, val_func = res
@unpack a_grid = par

# exercise 1.
# savings of the worker of 20 years old with high ans low productivity shocks
plot(a_grid, pol_func[:, 1, 20], label = "High productivity", xlabel = "Current assets, a", ylabel = "Future assets, a", title = "Savings of the worker at 20 years old")
plot!(a_grid, a_grid, label = "45 degree line", linestyle = :dash)
plot!(a_grid, pol_func[:, 2, 20], label = "Low productivity", legend =:bottomright)

savefig(figpath*"savings_at_20.png")

# value function of a retired of 50 years old
plot(a_grid, val_func[:, 1, 50], xlabel = "Current assets, a", ylabel = "Value", title = "Value function of the retired at the age 50")
savefig(figpath*"value_at_50.png")

# create table for exercise 3
using Latexify

benchmark_model = [1 2;
                   3 4;
                   5 6;
                   7 8;
                   9 10;
                   11 12;
                   13 14]
norisk_model = benchmark_model
exlab_model  = benchmark_model