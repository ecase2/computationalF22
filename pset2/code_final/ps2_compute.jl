using Parameters, Plots, Statistics #import the libraries we want
include("PS2_model.jl") #import the functions that solve our growth model

prim, res = Initialize(0.994) #initialize primitive and results structs

### Solve the model 
### Return net wealth
d = ClearMarket(prim, res)

### Compute the Gini coefficient
sum_wealth, sum_people = Lorenz(prim, res)
Gini_coef = Gini(prim, res)

# make graphs and lorenz curve
include("ps2_graphs.jl")


# @unpack val_func, pol_func, distr, q = res
# @unpack a_grid, s_grid = prim

# ############## Make plots
# ### distribution

# Plots.plot(a_grid .+ 1, res.distr[:, 1], title= "Distribution of wealth (a+s)", xlim = [-2, 3], label = "Good shock")
# Plots.plot!(a_grid .+ 0.5, res.distr[:, 2], title= "Distribution of wealth (a+s)", xlim = [-2, 3], xlabel = "Wealth", label = "Bad shock")
# Plots.savefig("PS2_Distribution.png")

# #value function
# Plots.plot(a_grid, val_func[:, 1], title="Value Function", label = "Good shock")
# Plots.plot!(a_grid, val_func[:, 2], title="Value Function", label = "Bad shock", yaxis = "Value", xaxis = "Current assets", legend =:bottomright)
# Plots.savefig("PS2_Value_Functions.png")

# #policy function
# Plots.plot(a_grid, pol_func[:, 1], title="Policy Functions", label = "Good shock")
# Plots.plot!(a_grid, pol_func[:, 2], title="Policy Functions",  label = "Bad shock", yaxis = "Future assets", xaxis = "Current assets", legend =:bottomright)
# plot!(a_grid, a_grid, label = "45* line", linestyle = :dash)
# Plots.savefig("PS2_Policy_Functions.png")

# #changes in policy function
# pol_func_δ_1 = copy(pol_func[:, 1]).-a_grid
# pol_func_δ_2 = copy(pol_func[:, 2]).-a_grid

# Plots.plot(a_grid, pol_func_δ_1, title="Policy Functions Changes", label = "Good shock")
# Plots.plot!(a_grid, pol_func_δ_2, title="Policy Functions Changes", label = "Bad shock")
# Plots.savefig("PS2_Policy_Functions_Changes.png")
