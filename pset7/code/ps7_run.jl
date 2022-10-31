#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 6
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

# import packages used to run the model
using Parameters, DataFrames

# import model functions
include("ps7_model.jl")


#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Task 1. Derivations - see write up.

### Case 1. Just identified case with 2 moments: mean and variance.
# Task 2. Simulate the true data.

# Need to run for each type

par, res = Initialize(; type = 2)
res.y = true_data(par, res)

# Task 3. Sumulate errors for simulations.
res.ε_grid = errors_sim(par, res)

# Task 4.
# (a) Plot and consistent estimate b1 with W = I
res.W = Diagonal(ones(par.n))

# Plot
x_ax, y_ax, z_ax = plot_J(par, res; ρl = 0.35, ρstep = 0.05, ρh = 0.65, σl = 0.8, σstep = 0.05, σh = 1.2)
Plots.plot(x_ax, y_ax, z_ax, st=:surface, camera=(-15, 15))

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

# (d) J-test WRONG
j_st = J_test(par)

# Task 5.
### Case 2. Just identified case with 2 moments: variance and autocovariance.


# Task 6.
### Case 3. Over-identified case with 3 moments: mean, variance and covariance.
