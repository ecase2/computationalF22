#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 5
    * AUTHORS:     Hanna Han, Emily Case, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps4_init.jl    initializes structs for all models
    *    ps4_model.jl   creates functions for models
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
include("HelpfulFunctions.jl")
include("ps5_model.jl")


#======================================================#
#       RUN MODELS
#------------------------------------------------------#

