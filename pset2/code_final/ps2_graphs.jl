using Printf 
figpath = joinpath(pwd(), "pset2/code_final/")
@unpack val_func, pol_func, distr, q = res
@unpack a_grid, s_grid = prim

qstar = @sprintf "%.6f" q
write(figpath*"qstar.tex", qstar)

# --------------------- II. --------------------
#value function
Plots.plot(a_grid, val_func[:, 1], title="Value Function", label = "Good shock")
Plots.plot!(a_grid, val_func[:, 2], title="Value Function", label = "Bad shock", yaxis = "Value", xaxis = "Current assets", legend =:bottomright)
Plots.savefig(figpath*"PS2_Value_Functions.png")

# (a) policy functions g(a,s) for s = e and s = u 
Plots.plot(a_grid, pol_func[:, 1], title="Policy Functions", label = "Good shock")
Plots.plot!(a_grid, pol_func[:, 2], label = "Bad shock", yaxis = "Future assets", xaxis = "Current assets", legend =:bottomright)
plot!(a_grid, a_grid, label = "45* line", linestyle = :dash)
Plots.savefig(figpath*"PS2_Policy_Functions.png")

#changes in policy function
pol_func_δ_1 = copy(pol_func[:, 1]).-a_grid
pol_func_δ_2 = copy(pol_func[:, 2]).-a_grid

Plots.plot(a_grid, pol_func_δ_1, title="Policy Functions Changes", label = "Good shock")
Plots.plot!(a_grid, pol_func_δ_2, title="Policy Functions Changes", label = "Bad shock")
Plots.savefig(figpath*"PS2_Policy_Functions_Changes.png")

# (b) cross-sectional distribution of wealth for those employed and those unemployed on the same graph.
Plots.plot(a_grid .+ 1, res.distr[:, 1], title= "Distribution of wealth (a+s)", xlim = [-2, 3], label = "Good shock")
Plots.plot!(a_grid .+ 0.5, res.distr[:, 2], title= "Distribution of wealth (a+s)", xlim = [-2, 3], xlabel = "Wealth", label = "Bad shock")
Plots.savefig(figpath*"PS2_Distribution.png")



# --------------------- III. --------------------

#(a) Plot λ(a,s)
λ = calcCE(prim, res)

Plots.plot(a_grid, λ[:,1], title = "Consumption equivalent λ", label = "Employed")
Plots.plot!(a_grid, λ[:,2], title = "Consumption equivalent", label = "Unemployed", legend =:bottomright)
Plots.savefig(figpath*"PS2_IIIa.png")

#(b) Find Welfare_FB, Welfare_Inc, and WG
W_FB = @sprintf "%.4f" calcWelfareFB(prim)
W_Inc = @sprintf "%.4f" calcWelfareInc(prim, res)
W_WG = @sprintf "%.4f" calcWelfareGain(prim, res)
write(figpath*"W_FB.tex", W_FB)
write(figpath*"W_Inc.tex", W_Inc)
write(figpath*"W_WG.tex", W_WG)

#(c) Find fraction of population that benefit from complete markets
frac = @sprintf "%.1f" voteforcompmarkets(res)*100
write(figpath*"popfrac_voting.tex", frac)