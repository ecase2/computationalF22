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
if pwd() == "C:\\Users\\79267"
    cd("C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22")
end

root     = joinpath(pwd(), "pset4")
codepath = joinpath(root, "code")
figpath  = joinpath(root, "figures")

# import packages used to run the model
using Parameters, DataFrames, CSV, Statistics

# import model functions
include("ps4_init.jl")
include("ps4_model.jl")
include("ps4_figures.jl")
include("ps3_model.jl")


#======================================================#
#       RUN MODELS
#------------------------------------------------------#

# find old and new steady state objects
SS_old = SolveModel()
SS_new = SolveModel(θ = 0.0)

# exercise 1
@time prim, res, K, L, EV, voters = solveModel(SS_old, SS_new)

# exercise 2
prim2, res2, K2, L2, EV2, voters2 = solveModel(SS_old, SS_new, t_noss = 21)


#======================================================#
#       WRITE AND SAVE GRAPHS
#------------------------------------------------------#

rPath, wPath, LPath, KPath = graphPath(prim, res)
rPath2, wPath2, LPath2, KPath2 = graphPath(prim2, res2, cf = "2")

EVgraph = graphEV(prim, EV)
EVgraph2 = graphEV(prim, EV2, cf = "2")
