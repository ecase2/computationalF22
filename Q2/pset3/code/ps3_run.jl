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
function read_data(; file_data::String, file_instruments::String, file_sims::String)

    # Car characteristics
    df = DataFrame(load(file_data))
    df = df .+ 0.0

    N = length(df[:, 1])      # number of observations
    Nx = length(df[1, :])     # number of car characteristics

    # Instruments
    iv = DataFrame(load(file_instruments))
    iv = iv .+ 0.0
    Niv = length(iv[1, :])      # number of observations
#    Nx = length(df[1, :])     # number of car characteristics

    # Simulations (random coefficients)
    sim = DataFrame(load(file_sims))
    sim = sim .+ 0.0
    Nsim = length(sim[1, :])

    return df, N, Nx, iv, Niv, sim, Nsim
end



# Read the data
df, N, Nx, iv, Niv, sim, Nsim = read_data(; file_data = file_car_char, file_instruments = file_iv, file_sims = file_sim)

######################## TRANSLATE JF's OX CODE "car_demand_ps" INTO JULIA ########################
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
mEta = Matrix(sim)
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


### Let save some objects in the results structure for simplicity - will use these objects in some functions later
mutable struct Results
    #aMu::Vector{Float64}
    vShare::Vector{Float64}
    aProductID::Vector{Vector{Int64}}
    aZ::Vector{Matrix{Float64}}
    mEta::Matrix{Float64}
    Sim::Int64
    Share::Float64
    T::Int64
end

res = Results(vShare, aProductID, aZ, mEta, Sim, 0, T)

function value(res::Results; vParam::Float64, t::Int64)
    @unpack aProductID, Sim, aZ  = res

    rowid = aProductID[t]

    mMu = zeros(length(rowid), Sim) # mu for each product and each shock
    for j = 1:length(rowid) # for each possible price (product)
        for i = 1:Sim # for each possible shock
            mMu[j, i] = vParam*aZ[t][j, i]
            mMu[j, i] = exp(mMu[j, i])
        end
    end
    return mMu
end

mMu = value(res; vParam = 0.6, t = 1)


function demand_Newton(res::Results; aDel::Vector{Float64}, t::Int64, vParam::Float64)
    @unpack aProductID, Share = res

    aShare = 0

    rowid = aProductID[t] # choose observations (products) for a given year

    mMu = value(res; vParam = vParam, t = t)

    Del = aDel[rowid]

    eV = exp.(Del).*mMu
    mS = eV./(1 .+ sum(eV, dims = 1))
    vShat = mean(mS)

    mD = zeros(length(rowid), length(rowid))
    for j = 1:length(rowid), k = 1:length(rowid)
        if j == k
            mD[j, k] = mean(mS[j, :].*(1 .- mS[j, :]))
        else
            mD[j, k] = mean(-mS[j, :].*mS[k, :])
        end
    end
    Jacobian = mD

    aShare = vShat

    return aShare, Jacobian
end


function demand_contraction(res::Results; aDel::Vector{Float64}, t::Int64, vParam::Float64)
    @unpack aProductID, Share = res

    aShare = 0

    rowid = aProductID[t]

    mMu = value(res; vParam = vParam, t = t)
    Del = aDel[rowid]

    eV = exp.(Del).*mMu
    mS = eV./(1 .+ sum(eV, dims = 1))
    vShat = mean(mS)

    aShare = vShat

    return aShare
end

function inverse(res::Results; aDel::Vector{Float64}, eps1::Float64, eps::Float64, vParam::Float64, maxit::Int64 = 1000)
    @unpack T, vShare = res

    vShat = 0.54

    vIT = 0

    maxT = T

    temDelta = aDel

    for t = 1:maxT
        println("t = $t")
        mu = value(res; vParam, t)
        rowid = aProductID[t]
        vIT =  0
        f = 1000
        err = maximum(abs.(f))

        while (err > eps && vIT < maxit)
            vIT += 1

            if (err > eps1)
                vShat = demand_contraction(res; aDel = temDelta, t = t, vParam = vParam)
                f = log.(vShare[rowid]) .- log(vShat)  #/* Zero function: Log difference in shares */
                if vIT < length(rowid)-1
                    temDelta[rowid] = temDelta[rowid] .+ f # /* Contraction mapping step */
                end
            elseif (err <= eps1)
                vShat, Jac = demand_Newton(res; aDel = temDelta, t = t, vParam = vParam)
                f = log.(vShare[rowid]) .- log.(vShat)# /* Zero function: Log difference in shares */
                temDelta[rowid] = temDelta[rowid] .+ inv(Jac/vShat)*f #/* Newton step */
                println("Jac = $Jac")
            end
            err = maximum(abs.(f))
        end
    end
    return temDelta
end


temD = inverse(res; aDel = vDelta0, eps1 = 1.0, eps = 10e-12, vParam = 0.6, maxit = 1000)

# Have to rewrite this function since it will be optimized?
function gmm_obj(res::Results; vDel::Vector{Float64}, const adFunc, const avScore, const amHessian)
    @unpacl vDelta0 = res
    # Invert demand
    eps1 = 1
    eps = 10^(-12)
    inverse(res; vDel = vDelta0, eps1 = eps1, eps = eps, iprint = iprint, vParam = x, maxit = 1000)

    # Quality decomposition */
    decl vLParam=ivreg(vDelta0,mX,mIV,A);
    decl vXi=vDelta0-mX*vLParam;
    # GMM objective function */
    mG=vXi.*mIV;
    decl scale=100;
    decl g=sumc(mG);
    if(isnan(vDelta0)==1) adFunc[0]=.NaN;
    else adFunc[0]=double(-g*A*g'/scale);
    return 1;
}



d = inverse(res; vDel = vDelta0, eps1 = 0.0, eps = 10e-12,  iprint = 0, vParam = 0.6, maxit = 100)
vDelta0 - d
vDelta0


mu = value(res)
vDelta0
sh, jac = demand(res; mMu = mu, aJac = [0.0 0.0; 0.0 0.0], vDelta = vDelta0, t = 1, vParam = 0.6)
