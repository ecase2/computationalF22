#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 7
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps6_model.jl        creates functions for models
=#

#======================================================#
#       SET UP
#------------------------------------------------------#

# define directory paths
#   NOTE: in visual studio code, make sure you have the computationalF22 folder opened, so that
#   pwd() automatically returns the file path to that folder
if pwd() == "C:\\Users\\79267"
    cd("C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22")
end

root     = joinpath(pwd(), "pset7")
codepath = joinpath(root, "code")
figpath = joinpath(root, "figures")

# import packages used to run the model
using Parameters, Distributions, Random, Optim, PyPlot, Plots, LinearAlgebra, StatsBase, Latexify

# import model functions
include("ps7_model.jl")
include("ps7_figures.jl")


#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Have to run the lines 44 - 103 three times:
# for type = 1 (use mean and variance)
# for type = 2 (use variance and autocorrelation)
# for type = 3 (use mean, variance, autocorrelation)

type = 1

if type == 1
    c = "type1"
    name = "type1.png"
elseif type == 2
    c = "type2"
    name = "type2.png"
elseif type == 3
    c = "type3"
    name = "type3.png"
end

par, res = Initialize(; type = type)
res.y, res.m_data, res.var_data, res.acor_data = true_data(par, res) # it is not dependent on the type (case) and will be the same for all types.

# Task 3. Sumulate errors for simulations.
res.ε_grid = errors_sim(par, res) # will be the same for all types

# Task 4.
# (a) Plot and consistent estimate b1 with W = I
res.W = Diagonal(ones(par.n))

# Plot
x_ax, y_ax, z_ax = plot_J(par, res; ρl = 0.35, ρstep = 0.05, ρh = 0.65, σl = 0.8, σstep = 0.05, σh = 1.2)
Plots.plot(x_ax, y_ax, z_ax, st=:surface, camera=(-15, 15))
Plots.savefig(joinpath(figpath, name))

# Estimates
sol = optimize(J, [0.5, 1.0])
res.b11, res.b21 = sol.minimizer[1], sol.minimizer[2]

# (b) Find W* and efficient estimate b2 with W = W*
dif = dif_for_gamma(par, res)
Γ0, Γj = gamma(par, res; dif = dif)
S_hat = S(par, res; Γ0 = Γ0, Γj = Γj)
res.W = inv(S_hat)

# Estimates
sol2 = optimize(J, [0.5, 1.0])
res.b12, res.b22 = sol2.minimizer[1], sol2.minimizer[2]

# (c) Derivatives and standard errors
∇ = derivatives(par, res)
errors = st_errors(par; ∇ = ∇, S_hat = S_hat)

# (d) J-test
j_st = J_test(par)

out = [res.b11, res.b21, res.b12, res.b22, res.W, ∇, errors, j_st]
estimates_table(; output1 = [out[1], out[2]], output2 = [out[3], out[4]], cf = c)
errors_table(; output1 = [out[3], out[4]], output2 = out[7], cf = c)
jacobian_table(; output1 = out[6], cf = c)

# save J-tests
if type == 1
    j1 = out[end]
elseif type == 2
    j2 = out[end]
elseif type == 3
    j3 = out[end]
end

# Bootstrap part
par, res = Initialize(; type = 3)
b1, b2 = bootstrap(par, res; N = 100000)

graph(; x = b2[:, 1], param = "rho")
graph(; x = b2[:, 2], param = "sigma")
