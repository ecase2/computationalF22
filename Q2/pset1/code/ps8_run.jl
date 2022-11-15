#=
    * PROJECT:     COMPUTATIONAL FALL 2022 PSET 8
    * AUTHORS:     Hanna Han, Anna Lukianova
    * CONTENTS:    THIS FILE RUNS THE ENTIRE PROGRAM.
    * OTHER FILES:
    *    ps6_model.jl        creates functions for models
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

root     = joinpath(pwd(), "pset8")
codepath = joinpath(root, "code")

# import packages used to run the model

# Pkg.add("StatFiles")
using Parameters, DataFrames, StatFiles, ForwardDiff

# import model functions
include("ps8_model.jl")

#======================================================#
#       RUN
#------------------------------------------------------#
# Run model

# Prepare the data
# Read the dataset
# Anna: I could not to put this function in the model code.
function prepare_data()
    # Anna: have to put the entire path...
    df = DataFrame(load("C:\\Users\\79267\\Documents\\UW_PhD\\3rd_year\\computationalF22\\pset8\\code\\Mortgage_performance_data.dta"))

    # Observed characteristics
    X = copy(df)
    X = select!(X,  ["i_large_loan", "i_medium_loan", "rate_spread", "i_refinance", "age_r",  "cltv",
                "dti", "cu", "first_mort_r", "score_0", "score_1", "i_FHA",
                 "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"])
    X = Matrix(X)
    # Add the column for constants
    X_add = ones(length(X[:, 1]))
    X = hcat(X, X_add)

    # Dependendent binary variable
    Y = copy(df)
    Y = select!(Y, "i_close_first_year")
    Y = Matrix(Y)

    N = size(Y)[1] # number of observations
    K =  size(X)[2] # number of independent variables

    return Y, X, N, K
end


Y, X, N, K = prepare_data()
res = Initialize(; N = N, K = K)

# Task 1. Evaluate log-likelihood function, score of log-likelihood function and Hession
β = vcat(-1, zeros(K-1))

LL = LogLikelihood(β)
scoreLL = score(β)
H = hessian(β)

# Task 2. Compare evaluated log-likelihood function, score of log-likelihood function and Hession
# with numerical ones
res = Initialize(; N = N, K = K)

NumScore = ForwardDiff.gradient(LogLikelihood, β)

dif = sqrt(dot(scoreLL - NumScore, scoreLL - NumScore))
# Anna: The difference is just 0.0 - can it be true?

# The Hessian
res = Initialize(; N = N, K = K)

NumH = ForwardDiff.hessian(LogLikelihood, β)

difH = sum(abs.(H .- NumH))

# Task 3. Solution of the maximum likelihood problem using a Newton algorithm.
@time MLL = NewtonMethod(β; K = K)
MLL

# Task 4. Compare the solution and numerical speed with two optimization packages provided
# by your favorite software: BFGS and Simplex.
guess_init = zeros(K)
guess_init[1] = 1
opt = optimize(LogLikelihood, guess_init)

optimize(LogLikelihood, guess_init, Optim.Options(iterations = 1000))

# Produce output table
write_table(output)
