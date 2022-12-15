#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 3 (Q2)
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

# The code was written based on JF's sample code

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

    # Simulations (random coefficients)
    sim = DataFrame(load(file_sims))
    sim = sim .+ 0.0
    Nsim = length(sim[1, :])

    return df, N, Nx, iv, Niv, sim, Nsim
end

# Results structure - will save object which we will use in some functions later
mutable struct Results
    #aMu::Vector{Float64}
    vShare::Vector{Float64}
    aProductID::Vector{Vector{Int64}}
    aZ::Vector{Matrix{Float64}}
    mEta::Matrix{Float64}
    Sim::Int64
    Share::Float64
    T::Int64
    vDelta0::Vector{Float64}
end


# IV-regression functions
function ivreg(; mY::Vector, mVariables::Matrix, mInstruments::Matrix, mWeight::Matrix)

    mInvXZX = inv((transpose(mVariables)*mInstruments)*mWeight*(transpose(mInstruments)*mVariables))

    if mInvXZX == 0
        vIV = ones(length(mVariables[1, :]), 1)*NaN
    else
        vIV = mInvXZX*(transpose(mVariables)*mInstruments)*mWeight*(transpose(mInstruments)*mY)
    end

    return vIV
end


# Pre-compute the row IDs for each market (each market = each year?)
function id_markets()
    add1 = [findall(x-> x == vYear[1], df[:, "Year"])]
    aProductID = Vector{Vector{Int64}}(add1)

    for i = 2:T
        add = findall(x-> x == vYear[i], df[:, "Year"])
        push!(aProductID, add)
    end

    return aProductID
end

# Pre-compute interaction between price and random-coefficient
# Two dimenional arrays of JxSim matrices (J is number of products in a given year and Sim is number of possible shocks):

# aZ include T matrices

function price_ran_coeff_inter(; aProductID::Vector{Vector{Int64}})
    aZ = Vector{Matrix{Float64}}([])

    for i = 1:T
        k = mZ[aProductID[i]]
        add = zeros(length(k), length(mEta[:, 1]))
        for j = 1:length(k)
            add[j, :] = k[j]*mEta[:, 1]
        end
        push!(aZ, add)
    end

    return aZ
end


# This function evaluates the idiosyntractic component of utility
# vParam is just one parameter λ

function value(res::Results; vParam::Float64, t::Int64, rowid::Vector{Int64})
    @unpack aProductID, Sim, aZ  = res

    mMu = zeros(length(rowid), Sim) # mu for each product and each shock
    for j = 1:length(rowid) # for each possible price (product)
        for i = 1:Sim # for each possible shock
            mMu[j, i] = vParam*aZ[t][j, i]
            mMu[j, i] = exp(mMu[j, i])
        end
    end
    return mMu
end

function demand_Newton(res::Results; aDel::Vector{Float64}, t::Int64, vParam::Float64, rowid::Vector{Int64})
    @unpack aProductID, Share = res

    aShare = 0

    mMu = value(res; vParam = vParam, t = t, rowid = rowid)

    Del = aDel

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


function demand_contraction(res::Results; aDel::Vector{Float64}, t::Int64, vParam::Float64, rowid::Vector{Int64})
    @unpack aProductID, Share = res

    aShare = 0

    mMu = value(res; vParam = vParam, t = t, rowid = rowid)
    Del = aDel

    eV = exp.(Del).*mMu
    mS = eV./(1 .+ sum(eV, dims = 1))

    vShat = mean(mS)

    aShare = vShat

    return aShare
end


function inverse(res::Results; eps1::Float64, eps::Float64, vParam::Float64, maxit::Int64 = 1000, t::Int64)
    @unpack aProductID, T, vShare, vDelta0 = res

    vShat = 0.54

    vIT = 0

    maxT = T

        rowid = aProductID[t]

        temDelta = vDelta0[rowid]
        println("initial temDelta[1:5] = ", temDelta[1:5])
        share = vShare[rowid]

        vIT =  0
        f = 100*ones(length(rowid))

        err = sqrt(sum(f.^2))

        while (err > eps && vIT < maxit)
            vIT += 1

            if (err > eps1)
                vShat = demand_contraction(res; aDel = temDelta, t = t, vParam = vParam, rowid = rowid)
                println("Contraction: it = $vIT, vShat = $vShat")
                f = log.(share) .- log(vShat)  #/* Zero function: Log difference in shares */
                temDelta = temDelta .+ f # /* Contraction mapping step */
            elseif (err <= eps1)
                vShat, Jac = demand_Newton(res; aDel = temDelta, t = t, vParam = vParam, rowid = rowid)
                f = log.(share) .- log.(vShat)# /* Zero function: Log difference in shares */
                temDelta = temDelta .+ inv(Jac/vShat)*f #/* Newton step */
            end
            err = sqrt(sum(f.^2))
            println("err = $err")
        end

    return temDelta
end
