# set some defaults for the plots:
default(titlefont = (20, "times"), linewidth = 2)

# Anna: Not sure whether I need to include some files here. Probably, not.

@unpack pol_func, val_func = res
@unpack a_grid = par

# exercise 1.
# savings of the worker of 20 years old with high ans low productivity shocks
# Anna: I don't like the policy function
StatsPlots.plot(a_grid, pol_func[:, 1, 20], label = "High productivity", xlabel = "Current assets, a", ylabel = "Future assets, a", title = "Savings of the worker at 20 years old")
StatsPlots.plot!(a_grid, pol_func[:, 2, 20], label = "Low productivity", legend =:bottomright)
savefig(joinpath(figpath, "pol_func.png"))

# value function of a retired of 50 years old
StatsPlots.plot(a_grid, val_func[:, 1, 50], xlabel = "Current assets, a", ylabel = "Value", legend = false, title = "Value function of the retired at the age 50")
savefig(joinpath(figpath, "savings_at_20.png"))

# value function of a retired of 50 years old
StatsPlots.plot(a_grid, val_func[:, 1, 50], xlabel = "Current assets, a", ylabel = "Value", title = "Value function of the retired at the age 50")
savefig(joinpath(figpath, "value_at_50.png"))

plot(a_grid, pol_func[:, 1, 20], label = "High productivity", xlabel = "Current assets, a", ylabel = "Future assets, a", title = "Savings of the worker at 20 years old")
plot!(a_grid, a_grid, label = "45 degree line", linestyle = :dash)
plot!(a_grid, pol_func[:, 2, 20], label = "Low productivity", legend =:bottomright)

savefig(joinpath(figpath, "savings_at_20.png"))

# value function of a retired of 50 years old
plot(a_grid, val_func[:, 1, 50], xlabel = "Current assets, a", ylabel = "Value", title = "Value function of the retired at the age 50")
savefig(joinpath(figpath,"value_at_50.png"))

# create table for exercise 3
using Latexify
set_default(fmt = "%.2f", convert_unicode = false)

x = hcat(bm_ss, bm_ss, bm_ss, bm_ss, bm_ss, bm_ss)

tabhead = [ "with SS", "wo SS", "with SS", "wo SS","with SS", "wo SS"]
tabside = [ "capital K", "labor L", "wage w", "interest r", "pension benefit b"]
temptable = latexify(x, env = :tabular, latex = false, side = tabside, head = tabhead)

temptable[1:26]
insert = "& \\multicolumn{2}{Benchmark model} & \\multicolumn{2}{No risk, \$z^L = z^H = 0.5\$} & \\multicolumn{2}{Exogenous labor, \$ \\gamma = 1\$} \n"
write(joinpath(figpath, "resultstable.tex"), temptable[1:26]*insert*temptable[27:end])