using Distributed
@everywhere using Parameters, Plots, SharedArrays #import the libraries we want
@everywhere include("model.jl") #import the functions that solve our growth model


prim, res = Initialize() #initialize primitive and results structs
runtime = @elapsed Solve_model(prim, res) #solve the model and save time 

addprocs(5) 
runtime_par = @elapsed Solve_model(prim, res) #solve the model and save time 


@unpack val_func, pol_func = res
@unpack k_grid = prim

##############Make plots and tex snippets 
figpath = joinpath(pwd(), "pset1/figures")
write(joinpath(figpath, "julia_runtime.tex"), string(runtime))
write(joinpath(figpath, "julia_runtime_parallel.tex"), string(runtime_par))
#value function
plot(k_grid, val_func[:,1], title="Value Functions", xlabel = "K", ylabel = "V(K,Z)", label = "Zg")
plot!(k_grid, val_func[:,2], label = "Zb")
savefig(figpath*"/valfunctions.png")

#policy functions
plot(k_grid, pol_func[:,1], title="Policy Functions", xlabel = "K", ylabel = "K'", label = "Zg")
plot!(k_grid, pol_func[:,2], label = "Zb")
plot!(k_grid, k_grid, label = "45 degree line")
savefig(figpath*"/polfunctions.png")

#changes in policy function
pol_func_δ_g = copy(pol_func[:,1]).-k_grid
pol_func_δ_b = copy(pol_func[:,2]).-k_grid
plot(k_grid, pol_func_δ_g, title="Policy Functions Changes", xlabel = "K", label = "Zg")
plot!(k_grid, pol_func_δ_b, label = "Zb")
savefig(figpath*"polfuncchanges.png")