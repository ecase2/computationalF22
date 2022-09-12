using Parameters, Plots #import the libraries we want
include("PS1_model.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func = res
@unpack k_grid, z_grid = prim

##############Make plots
#value function
Plots.plot(k_grid, val_func[:, 1], title="Value Function", label = "Good shock")
Plots.plot!(k_grid, val_func[:, 2], title="Value Function", label = "Bad shock", legend =:bottomright)
Plots.savefig("PS1_Value_Functions.png")

#policy functions
Plots.plot(k_grid, pol_func[:, 1], title="Policy Functions", label = "Good shock")
Plots.plot!(k_grid, pol_func[:, 2], title="Policy Functions",  label = "Bad shock", legend =:bottomright)
Plots.savefig("PS1_Policy_Functions.png")

#changes in policy function
pol_func_δ_1 = copy(pol_func[:, 1]).-k_grid
pol_func_δ_2 = copy(pol_func[:, 2]).-k_grid

Plots.plot(k_grid, pol_func_δ_1, title="Policy Functions Changes. Good shock")
Plots.savefig("PS1_Policy_Functions_Changes_good.png")

Plots.plot(k_grid, pol_func_δ_2, title="Policy Functions Changes. Bad shock")
Plots.savefig("PS1_Policy_Functions_Changes_bad.png")
