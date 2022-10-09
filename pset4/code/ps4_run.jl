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
#include("ps4_model.jl")
include("ps4_model_grids.jl")
include("ps4_figures.jl")
include("ps3_model.jl")
### Need to think how to add old code
#======================================================#


#======================================================#
#       RUN MODELS
#------------------------------------------------------#

### find old and new steady state objects
res_old = SolveModel(; al = 0.0, au = 20.0, na = 200, θ = 0.11, z = [3.0, 0.5], γ = 0.42, T = 1, t_noss = 2, λ = 0.5, iter = 1000, tol = 0.005)
res_new = SolveModel(; al = 0.0, au = 20.0, na = 200, θ = 0.0, z = [3.0, 0.5], γ = 0.42, T = 1, t_noss = 2, λ = 0.5, iter = 1000, tol = 0.005)

res_old.Γ

# exercise 1
# find transition path
@time prim, res, K, L = solveModel(res_old, res_new; al = 0.0, au = 20.0, na = 200, θ = 0.11, z = [3.0, 0.5], γ = 0.42, T = 50, t_noss = 2, tol = 0.1, iter = 1000, λ = 0.90)

# exercise 2
prim2, res2, K2, L2  = solveModel(res_old, res_new; al = 0.0, au = 20.0, na = 200, θ = 0.11, z = [3.0, 0.5], γ = 0.42, T = 50, t_noss = 21, tol = 0.1, iter = 1000, λ = 0.90)

#======================================================#


#======================================================#
#       WRITE AND SAVE GRAPHS
#------------------------------------------------------#

rPath, wPath, LPath, KPath = graphPath(res; T = 50)
rPath2, wPath2, LPath2, KPath2 = graphPath(res2, ins2, cf = "2")

EVgraph2 = graphEV(prim, ev)
EVgraph2 = graphEV(prim, ev2; cf = "")

#======================================================#
