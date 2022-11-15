#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 2 (Q2)
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

mutable struct Results
    Λ::Array{Float64, 1}
    LL::Float64
    scoreLL::Array{Float64,1}
    Hessian::Array{Float64,2}
    β::Array{Float64, 1}
end

### Put initialization somewhere else...
function Initialize(; N::Int64, K::Int64)
    Λ = zeros(N)
    LL = 0.0
    scoreLL = zeros(K)
    Hessian =  zeros(K, K)
    β = ones(K).*10
    res = Results(Λ, LL, scoreLL, Hessian, β)
    return res
end


# Exercise 1
function LogLikelihood(β)
    @unpack LL = res
#    @unpack N, K = par
    for i in 1:N
        LL = LL .+ log((exp(dot(X[i, :], β))/(1 + exp(dot(X[i, :], β))))^Y[i,1]*(1-(exp(dot(X[i, :], β))/(1 + exp(dot(X[i, :], β)))))^(1-Y[i,1]))
    end
    return LL
end

function score(β::Vector{Float64})
    @unpack Λ = res

    for i = 1:N
        res.Λ[i] = exp(X[i, :]'*β)/(1+exp(X[i, :]'*β))
        res.scoreLL += (Y[i] - Λ[i])*X[i, :]
    end

    return res.scoreLL

end

function hessian(β::Vector{Float64})
    @unpack Λ = res

    for i = 1:N
        res.Hessian += -Λ[i]*(1-Λ[i])*X[i, :]*transpose(X[i, :])
    end

    return res.Hessian

end


# Exercise 2. Compare the score and hessian obtained from (1) with the numerical first and second
# derivative of the log-likelihood
# See in the ps8_run.jl

# Exercise 3. Routine that solves the maximum likelihood problem using a Newton algorithm.
########### Exercise 3. A Newton Algorithm to solve the maximum lilelihood problem


function NewtonMethod(β, tol = 10e-2; K::Int64)
    @unpack Hessian, scoreLL = res

    diff = 10.0
    iter = 1

    β1 = zeros(K)

    while diff > tol && iter < 1000
        scoreLL = score(β)
        Hessian = hessian(β)
    #    println("β  = $β\n")
        β1 =  β - inv(Hessian)*scoreLL
        diff = sqrt(dot(β1 - β, β1 - β))
        if diff > tol
            β = β1
        #    println("Iter = $iter, diff = ", diff)
        end
        if diff <= tol
        #   println("Solution is found. That is β = ", β)
        end
        iter = iter + 1
    end
    return β1
end


# Exercise 4. Compare
# see ps8_run
