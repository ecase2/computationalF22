#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 3
    * AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps3_initialize.jl    initializes structs for all models
    *    ps3_model.jl         contains the main functions to run the model
    *    ps3_figs.jl          creates all graphs, figures, and .tex inputs
=#


#======================================================#
#       SET UP
#------------------------------------------------------#

# define directory paths
root     = joinpath(pwd(), "pset3")
codepath = joinpath(root, "code")
figpath  = joinpath(root, "figures")

# import packages used to run the model (plots and latexify are included in ps3_figs.jl)
using Parameters, DataFrames, CSV, Statistics, Plots, Latexify

# import model functions
include("ps3_initialize.jl")
include("ps3_model.jl")
include("ps3_figs.jl")
#======================================================#

#======================================================#
#       RUN MODELS
#------------------------------------------------------#

# Run benchmark model and create graphs
bm_ss, par, res   = SolveModel()
createAllGraphs(par, res)
bm_noss = SolveModel(θ = 0.0)

# No productivity shocks - this one has trouble
noshock_ss   = SolveModel(z = [0.5, 0.5], λ = 0.3)
noshock_noss = SolveModel(θ = 0.0, z = [0.5, 0.5], λ = 0.1)

# Inelastc labor supply
inelasticl_ss   = SolveModel(γ = 1.0)
inelasticl_noss = SolveModel(θ = 0.0, γ = 1.0)
#======================================================#

#======================================================#
#       WRITE TABLES AND GRAPHS
#------------------------------------------------------#
# NOTE: these need to be manually uploaded to the pset3/figures folder on Overleaf
resultsTable(bm_ss, bm_noss, noshock_ss, noshock_noss, inelasticl_ss, inelasticl_noss)
#======================================================#
