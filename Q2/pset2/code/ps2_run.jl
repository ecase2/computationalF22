#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 2 (Q2)
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps2_model.jl        creates functions for models
=#

#======================================================#
#       SET UP
#------------------------------------------------------#
# define directory paths
if pwd() == "C:\\Users\\79267"
    cd("C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22")
end

root     = joinpath(pwd(), "pset8")
codepath = joinpath(root, "code")

# import packages used to run the model
using Parameters, DataFrames, StatFiles, ForwardDiff

# import model functions
include("ps2_model.jl")


#======================================================#
#       RUN
#------------------------------------------------------#
