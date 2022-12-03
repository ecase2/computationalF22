#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 3 (Q2)
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

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



# IV-regression functions
function ivreg(; mY::Vector, mVariables::Matrix, mInstruments::Matrix, mWeight::Matrix)
    ### NOT SURE ABOUT SIZES OF THE FUNCTION ARGUMENTS

    mInvXZX = inv((transpose(mVariables)*mInstruments)*mWeight*(transpose(mInstruments)*mVariables))

    if mInvXZX == 0
        vIV = ones(length(mVariables[1, :]), 1)*NaN
    else
        vIV = mInvXZX*(transpose(mVariables)*mInstruments)*mWeight*(transpose(mInstruments)*mY)
    end

    return vIV
end



# This function evaluates the idiosyntractic component of utility

# aZ, aProductID should be in structure
# vParam is just one parameter λ
function value(res::Results; vParam::Float64 = 0.6, t::Int64 = 1)
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


function demand(res::Results; mMu::Matrix, aJac::Matrix{Float64}, vDelta::Vector{Float64}, t::Int64, vParam::Float64)
    @unpack aProductID, Share = res

    aShare = 0

    rowid = aProductID[t] # choose observations (products) for a given year

    eV = exp.(vDelta[rowid]).*mMu
    mS = eV./(1 .+ sum(eV, dims = 1)) ### NEED TO SUM ELEMENTS FROM COLUMNS (that is for fixed shock!) - choice probability of a given type
    vShat = mean(mS) # It is σ in our problem set - mean across all simulations

    if mean(aJac) != 0 # mD = Dσ (if we will use the Newton algorithm, need to compute the Jacobian)
        mD = zeros(length(rowid), length(rowid))
        for j = 1:length(rowid), k = 1:length(rowid)
            if j == k
                mD[j, k] = mean(mS[j, :].*(1 .- mS[j, :]))
            else
                mD[j, k] = mean(-mS[j, :].*mS[k, :])
            end
        end
        aJac = mD
    end
    #    mD = diag(meanr(mS.*(1-mS)))-setdiagonal(mS*mS'/Sim,0) ### NEED TO FIX THIS LINE

    aShare = vShat ### NEED TO BE CAREFUL

    return aShare, aJac # need to update res.Share with this value
end
