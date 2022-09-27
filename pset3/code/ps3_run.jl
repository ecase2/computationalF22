
#=
    COMPUTATIONAL FALL 2022 PSET 3
    AUTHORS: Hanna Han, Emily Case, Anna Lukianova
    THIS FILE RUNS THE ENTIRE PROGRAM. OTHER CODE FILES ARE:
    - ps3_initialize.jl: initializes structs for all models 
    - ps3_model.jl: contains the main functions to run the model
    - ps3_figs.jl: creates all graphs, figures, and .tex inputs
    NOTE TO HANNA AND ANNA: maybe we could have multiple model function files, that way each file is a little bit more concise. we can always include them one after the other in this file. it may also help with everyone working at the same time - emily 9.21.22
=#

# define directory paths
root     = joinpath(pwd(), "pset3")
codepath = joinpath(root, "code")
figpath  = joinpath(root, "code/figs/")

# import packages
using Parameters, DataFrames, CSV, Plots

# run model functions
include("ps3_initialize.jl")
include("ps3_model.jl")

#----------------------------------------------#
# run the benchmark model and two counterfactuals 
solve(modeltype = "benchmark")
solve(modeltype = "norisk")
solve(modeltype = "benchmark")



#----------------------------------------------#

# create graphs and other inputs for the writeup
# NOTE: these need to be manually uploaded to the pset3/figures folder on Overleaf
include("ps3_figs.jl")