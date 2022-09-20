####################
# 0. Set up
####################
using Parameters, Plots, Statistics #import the libraries we want
include("ps2_model.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func, distr, q = res
@unpack a_grid, s_grid = prim


####################
# 1. Question II
####################
I = Indicator(prim, res)
mean(I)
I[10, 11, :]
findall(x -> x == 1, I)

res.distr = Distr(prim, res)
sum(distr)
Plots.plot(a_grid, distr[:, 1], title= "Distribution of assets", label = "Good shock")
Plots.plot!(a_grid, distr[:, 2], title= "Distribution of assets", label = "Bad shock")

ED = ExcessDemand(prim, res)

d = ClearMarket(prim, res)

##############Make plots
#value function
Plots.plot(a_grid, val_func[:, 1], title="Value Function", label = "Good shock")
Plots.plot!(a_grid, val_func[:, 2], title="Value Function", label = "Bad shock", legend =:bottomright)
Plots.savefig("PS2_Value_Functions.png")

#policy functions
Plots.plot(a_grid, pol_func[:, 1], title="Policy Functions", label = "Good shock")
Plots.plot!(a_grid, pol_func[:, 2], title="Policy Functions",  label = "Bad shock", legend =:bottomright)
Plots.savefig("PS2 _Policy_Functions.png")

#changes in policy function
pol_func_δ_1 = copy(pol_func[:, 1]).-a_grid
pol_func_δ_2 = copy(pol_func[:, 2]).-a_grid

Plots.plot(a_grid, pol_func_δ_1, title="Policy Functions Changes. Good shock")
Plots.savefig("PS2_Policy_Functions_Changes_good.png")

Plots.plot(a_grid, pol_func_δ_2, title="Policy Functions Changes. Bad shock")
Plots.savefig("PS2_Policy_Functions_Changes_bad.png")


####################
# 2. Question III
####################
#(a) Plot λ(a,s)
λ = calcCE(prim, res)

Plots.plot(a_grid, λ[:,1], title = "Consumption equivalent", label = "Employed")
Plots.plot!(a_grid, λ[:,2], title = "Consumption equivalent", label = "Unemployed", legend =:bottomright)
Plots.savefig("PS2_IIIa.png")

#(b) Find Welfare_FB, Welfare_Inc, and WG
W_FB = calcWelfareFB(prim)
W_Inc = calcWelfareInc(prim, res)
WG = calcWG(prim, res)

#(c) Find fraction of population that benefit from complete markets
