#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 3 (Q2)
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps3_model.jl        creates functions for models
    *    ps3_figures.jl      creates tables
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

root     = joinpath(pwd(), "pset3")
codepath = joinpath(root, "code")
tabpath  = joinpath(root, "tables")

# using Parameters, Distributions, Random, Optim, PyPlot, Plots, LinearAlgebra, StatFiles, ForwardDiff, DataFrames, CSV

file_car_char = joinpath(codepath, "Car_demand_characteristics_spec1.dta")
file_iv = joinpath(codepath, "Car_demand_iv_spec1.dta")
file_sim = joinpath(codepath, "Simulated_type_distribution.dta")


# import packages used to run the model
using Parameters, DataFrames, StatFiles, ForwardDiff, Latexify, CSV, Distributions, Random, Optim

# import model functions
include("ps3_model.jl")
include("ps3_tables.jl")


#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Read the data
df, N, Nx, iv, Niv, sim, Nsim = read_data(; file_data = file_car_char, file_instruments = file_iv, file_sims = file_sim)

# Panel structure
n = N
vYear = sort(unique(df[:, "Year"]))
T = length(vYear)

# Outcome variables
vShare = df[:, "share"]
vDelta_iia = df[:, "delta_iia"]
vDelta0 = vDelta_iia
vPrice = df[:, "price"]

# Variables (all - exogenous and endogeneous)
varlist = ["price","dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990",
    "Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001",
    "Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012",
    "model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8","drive_2"]

# Only exogenous variables (price is excluded)
exo_varlist = ["dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990",
    "Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001",
    "Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012",
    "model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8"]

# All variables
mX = df[:, varlist]

# Load price and differentation IVs
aIVlist = ["i_import","diffiv_local_0","diffiv_local_1","diffiv_local_2","diffiv_local_3","diffiv_ed_0"]
mExclIV = iv[:, aIVlist]

# I AM NOT SURE ABOUT THE NEXT !!!!
#mIV = mPanelCharact[][find(aCharactName,exo_varlist)]~mExclIV;
mIV = hcat(df[:, exo_varlist], mExclIV)

# Non-linear attributes
mZ = df[:, "price"]


# Pre-compute the row IDs for each market (each market = each year?)
add1 = [findall(x-> x == vYear[1], df[:, "Year"])]
aProductID = Vector{Vector{Int64}}(add1)

for i = 2:T
    add = findall(x-> x == vYear[i], df[:, "Year"])
    push!(aProductID, add)
end

# Random coefficients
mEta = sim
Sim = length(mEta[:, 1])


# Pre-compute interaction between price and random-coefficient
# Two dimenional arrays of JxSim matrices (J is number of products in a given year and Sim is number of possible shocks):
# T x Nb of variables (what is Nb?)

# aZ include T matrices

aZ = Vector{Matrix{Float64}}([])

for i = 1:T
    k = mZ[aProductID[i]]
    add = zeros(length(k), length(mEta[:, 1]))
    for j = 1:length(k)
        add[j, :] = k[j]*mEta[:, 1]
    end
    push!(aZ, add)
end
