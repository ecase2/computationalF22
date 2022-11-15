#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 9
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps9_model.jl        creates functions for models
    *    ps9_tables.jl       creates functions for tables
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

root     = joinpath(pwd(), "pset9")
codepath = joinpath(root, "code")
tabpath  = joinpath(root, "tables")
# import packages used to run the model

# Pkg.add("StatFiles")
using Parameters, DataFrames, StatFiles, ForwardDiff

# import model functions
include("ps9_model.jl")
include("ps9_tables.jl")

#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Prepare the data
# Read the dataset from
if pwd() == "C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22"
    filename = "C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22\\pset9\\code\\Mortgage_performance_data.dta"
end

function prepare_data()
    df = DataFrame(load(filename))

    # Observed characteristics - time-invariant
    X = copy(df)
    X = select!(X,  ["score_0", "rate_spread", "i_large_loan", "i_medium_loan", "i_refinance", "age_r",  "cltv",
                "dti", "cu", "first_mort_r", "i_FHA",
                 "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"])
    X = Matrix(X)
    # Add the column for constants
    X_add = ones(length(X[:, 1]))
    X = hcat(X, X_add)

    # Observed characteristics - time-invariant
    Z = copy(df)
    Z = select!(Z,  ["score_0", "score_1", "score_2"])
    Z = Matrix(Z)

    # use Y = ["i_close_0", "i_close_1", "i_close_2"] to create a variable T
    Y = copy(df)
    Y = select!(I, ["i_close_0", "i_close_1", "i_close_2"])

    # construct T variable
    N = length(Y[:, 1]) # number of observations
    Kx = length(X[1, :])  # number of time-invariant observations
    Kz = length(Z[1, :]) # number of time-varying observations

    T = zeros(length(Y[:, 1]))

    for n = 1:N
        if Y[n, 1] == 1
            T[n] = 1
        elseif Y[n, 1] == 0 && Y[n, 2] == 1
            T[n] = 2
        elseif Y[n, 1] == 0 && Y[n, 2] == 0 && Y[n, 3] == 1
            T[n] = 3
        elseif Y[n, 1] == 0 && Y[n, 2] == 0 && Y[n, 3] == 0
            T[n] = 4
        end
    end
    return X, Z, T, Y, N, Kx, Kz
end


 X, Z, T, Y, N, Kx, Kz = prepare_data()
# res = Initialize(; N = N, K = K)
