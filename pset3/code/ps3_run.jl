
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
using Parameters, DataFrames, CSV, StatsPlots

# run model functions
include("ps3_initialize.jl")
include("ps3_model.jl")

#----------------------------------------------#
### RUN THE MODEL HERE USING MODEL FUNCTIONS ###
par, res, grid = Initialize()

fill_end_grids(par, res, grid) # it is better to put in somewhere else
backward_iteration(par, res, grid)

# to find distribution
F = distribution(par, res, grid)
sum(F) # it is less than one - bad...

# to find distribution using iterating procedure - works wrong too.
# NOTE: need to understand whether we use iteration to find the distribution or not (since now the model is finite)
# F_v2 = distribution_v2(par, res, grid)
# sum(F_v2)
# maximum(F_v2)


#----------------------------------------------#

# create graphs and other inputs for the writeup
# NOTE: these need to be manually uploaded to the pset3/figures folder on Overleaf
include("ps3_figs.jl")
