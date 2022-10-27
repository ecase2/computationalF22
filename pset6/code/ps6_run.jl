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
using Parameters, DataFrames

# import model functions
include("ps6_model.jl")


#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Task 1.
# Baseline
@time par1, res1 = solve_model(; cf = 10, α = 0)
#res1.p
#res1.pol_func

# With shocks. α = 1
@time par2, res2 = solve_model(; cf = 10, α = 1)
#res2.p
#res2.pol_func

# With shocks. α = 2
@time par3, res3 = solve_model(; cf = 10, α = 2)
#res3.p
#res3.pol_func


# Task 2


# Produce output table
write_table(output)
