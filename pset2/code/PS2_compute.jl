using Parameters, Plots, Statistics #import the libraries we want
include("PS2_model.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func, distr, q = res
@unpack a_grid, s_grid = prim

I = Indicator(prim, res)
mean(I)
I[10, 11, :]
findall(x -> x == 1, I)

res.distr = Distr(prim, res)
sum(res.distr)
Plots.plot(a_grid, res.distr[:, 1], title= "Distribution of assets", label = "Good shock")
Plots.plot!(a_grid, res.distr[:, 2], title= "Distribution of assets", label = "Bad shock")

ED = ExcessDemand(prim, res)

prim, res = Initialize(0.994) #initialize primitive and results structs
V_iterate(prim, res)
Distr(prim, res)
ExcessDemand(prim, res)
d = ClearMarket(prim, res)

for q = 0.7:0.05:0.9
    prim, res = Initialize(q)
    V_iterate(prim, res)
    println("pol_func[1:10, 1] = ", res.pol_func[1:10, 1], "q = $q")
    println("pol_func[1:10, 2] = ", res.pol_func[1:10, 2], "q = $q")
end 



##############Make plots
#value function
Plots.plot(a_grid, val_func[:, 1], title="Value Function", label = "Good shock")
Plots.plot!(a_grid, val_func[:, 2], title="Value Function", label = "Bad shock", legend =:bottomright) 

#policy function 
Plots.plot(a_grid, pol_func[:, 1], title="Policy Functions", label = "Good shock")
Plots.plot!(a_grid, pol_func[:, 2], title="Policy Functions",  label = "Bad shock", legend =:bottomright)
plot!(a_grid, a_grid, label = "45* line", linestyle = :dash)
Plots.savefig("PS2 _Policy_Functions.png")

#changes in policy function
pol_func_δ_1 = copy(pol_func[:, 1]).-a_grid
pol_func_δ_2 = copy(pol_func[:, 2]).-a_grid

Plots.plot(a_grid, pol_func_δ_1, title="Policy Functions Changes. Good shock")
Plots.savefig("PS2_Policy_Functions_Changes_good.png")

Plots.plot(a_grid, pol_func_δ_2, title="Policy Functions Changes. Bad shock")
Plots.savefig("PS2_Policy_Functions_Changes_bad.png")
