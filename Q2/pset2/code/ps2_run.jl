#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 2 (Q2)
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps2_model.jl        creates functions for models
    *    ps2_figures.jl      creates tables
=#


#======================================================#
#       SET UP
#------------------------------------------------------#

# define directory paths
#   NOTE: in visual studio code, make sure you have the computationalF22 folder opened, so that
#   pwd() automatically returns the file path to that folder
if pwd() == "C:\\Users\\79267"
    cd("C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22\\Q2")
end

root     = joinpath(pwd(), "pset2")
codepath = joinpath(root, "code")
tabpath  = joinpath(root, "tables")

# using Parameters, Distributions, Random, Optim, PyPlot, Plots, LinearAlgebra, StatFiles, ForwardDiff, DataFrames, CSV

filename = joinpath(codepath, "Mortgage_performance_data.dta")
file_weigths_dim1 = joinpath(codepath, "nodes_weights_dim1.csv")
file_weigths_dim2 = joinpath(codepath, "nodes_weights_dim2.csv")

# import packages used to run the model

# Pkg.add("StatFiles")
using Parameters, DataFrames, StatFiles, ForwardDiff, Latexify, CSV, Distributions

# import model functions
include("ps2_model.jl")
include("ps2_tables.jl")

#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Read the data
X, Z, T, Y, N, Kx, Kz = prepare_data(; filename)

par = Initialize(; Kx = Kx, Kz = Kz)

# Read nodes
nod1, nod21, nod22, w1, w2 = read_nodes(; file_weigths_dim1, file_weigths_dim2)


LL_quadr = quadr_LL(par; N = N, nodes1 = nod1, nodes21 = nod21, nodes22 = nod22,
    w1 = w1, w2 = w2, X = X, Z = Z, T = T)
