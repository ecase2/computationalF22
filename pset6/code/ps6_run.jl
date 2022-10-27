#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 6
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps6_model.jl        creates functions for models
    *    ps6_figures.jl      creates tables and figures
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

root     = joinpath(pwd(), "pset6")
codepath = joinpath(root, "code")
figpath  = joinpath(root, "figures")

# import packages used to run the model
using Parameters, DataFrames, CSV, Statistics, Plots, Latexify

# import model functions
include("ps6_model.jl")


#======================================================#
#       RUN
#------------------------------------------------------#
