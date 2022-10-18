#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 5
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps5_model.jl        creates functions for models
    *    HelpfulFunctions.jl contains helpful functions (written by Philip Coyle)
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

root     = joinpath(pwd(), "pset5")
codepath = joinpath(root, "code")

# import packages used to run the model
using Parameters, DataFrames, CSV, Statistics, Interpolations, Random, Distributions, Optim

# import model functions
include("HelpfulFunctions.jl")
include("ps5_model.jl")


#======================================================#
#       RUN MODELS
#------------------------------------------------------#

