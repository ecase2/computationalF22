#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 4
    * AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps4_init.jl    initializes structs for all models
    *    ps4_model.jl   creates functions for models
    *    ps4_figs.jl    creates functions for all graphs, figures, and .tex inputs
=#


#======================================================#
#       SET UP
#------------------------------------------------------#

# define directory paths
#   NOTE: in visual studio code, make sure you have the computationalF22 folder opened, so that
#   pwd() automatically returns the file path to that folder
root     = joinpath(pwd(), "pset4")
codepath = joinpath(root, "code")
figpath  = joinpath(root, "figures")

# import packages used to run the model
using Parameters, DataFrames, CSV, Statistics

# import model functions
include("ps4_init.jl")
include("ps4_model.jl")
include("ps4_figures.jl")
#======================================================#


#======================================================#
#       RUN MODELS
#------------------------------------------------------#

# exercise 1
ins, res = solveModel()

# exercise 2
ins2, res2 = solveModel()

#======================================================#


#======================================================#
#       WRITE AND SAVE GRAPHS
#------------------------------------------------------#

rPath, wPath, LPath, KPath = graphPath(res, ins)
rPath2, wPath2, LPath2, KPath2 = graphPath(res2, ins2)

#======================================================#
