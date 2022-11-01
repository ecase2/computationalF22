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

root     = joinpath(pwd(), "pset6")
codepath = joinpath(root, "code")
figpath = joinpath(root, "figures")

# import packages used to run the model
using Parameters, DataFrames, Plots, Latexify

# import model functions
include("ps6_model.jl")
include("ps6_figures.jl")


#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Task 1
@time par1, res1, output1 = solve_model(; cf = 10, α = 0) # Baseline
@time par2, res2, output2 = solve_model(; cf = 10, α = 1) # With shocks. α = 1
@time par3, res3, output3 = solve_model(; cf = 10, α = 2) # With shocks. α = 2

# Task 2
# Produce output table
write_table(output1, output2, output3)

# Task 3
# Produce plot
pol_func1 = res1.pol_func
pol_func2 = res2.pol_func
pol_func3 = res3.pol_func

graphExit(par1, pol_func1, pol_func2, pol_func3; plotname = "baseline")

# Task 4
# Run counterfactual
@time par1, res1, output12 = solve_model(; cf = 15, α = 0) # Baseline
@time par2, res2, output22 = solve_model(; cf = 15, α = 1) # With shocks. α = 1
@time par3, res3, output32 = solve_model(; cf = 15, α = 2) # With shocks. α = 2

pol_func12 = res1.pol_func
pol_func22 = res2.pol_func
pol_func32 = res3.pol_func

graphExit(par1, pol_func12, pol_func22, pol_func32; plotname = "counterfactual")

# If we show the results without α = 2
graphExit_no2(par1, pol_func12, pol_func22; plotname = "counterfactual")
