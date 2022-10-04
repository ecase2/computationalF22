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
#======================================================#


#======================================================#
#       RUN MODELS
#------------------------------------------------------#

prim, ins, res = initialize(0.11,  70, 0,  40,  200, 1)
# exercise 1
@time ins, res = solveModel(prim, ins, res; tol = 0.01, iter = 1000, λ = 0.995)

q1 = get_ss_distr(prim, ins, res)
sum(q1)
# exercise 2
ins2, res2 = solveModel()

#======================================================#


#======================================================#
#       WRITE AND SAVE GRAPHS
#------------------------------------------------------#

rPath, wPath, LPath, KPath = graphPath(res, ins)
rPath2, wPath2, LPath2, KPath2 = graphPath(res2, ins2, cf = "2")

EVgraph2 = graphEV(prim, ev)
EVgraph2 = graphEV(prim, ev2; cf = "")

#======================================================#
