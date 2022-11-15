#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 9
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

if pwd() == "C:\\Users\\79267"
    cd("C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22")
end

root     = joinpath(pwd(), "pset8")
codepath = joinpath(root, "code")


using Parameters, Distributions, Random, Optim, PyPlot, Plots, LinearAlgebra, StatFiles, ForwardDiff


mutable struct Results
end

### Put initialization somewhere else...
function Initialize(; N::Int64, Kx::Int64, Kz::Int64)

    res = Results( )
    return res
end
