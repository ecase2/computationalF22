# ECON899. Problem set 1.
# The code solves the neoclassical growth model with stochastic productivity.
# It is a modified version of the script "growth_julia" from Dean Corbae's website.
# Modifications are introduced by Emily Case (09/13/2022).
# This version allows for parallelization.


@everywhere @with_kw struct Primitives
    β::Float64 = 0.99                   #discount rate
    δ::Float64 = 0.025                  #depreciation rate
    α::Float64 = 0.36                   #capital share
    z::Vector{Float64} = [1.25, 0.2]    #zᵍ and zᵇ
    π::Matrix{Float64} = [[0.977, 0.074] [0.023, 0.926]]  #transition matrix
    k_min::Float64 = 0.01               #capital lower bound
    k_max::Float64 = 75.0               #capital upper bound
    nk::Int64 = 1000                    #number of capital grid points
    k_grid::Array{Float64,1} = collect(range(k_min, length = nk, stop = k_max)) #capital grid
end

#structure that holds model results
@everywhere mutable struct Results
    val_func::SharedArray{Float64, 2} #value function
    pol_func::SharedArray{Float64, 2} #policy function
end

#function for initializing model primitives and results
@everywhere function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = SharedArray{Float64}(zeros(prim.nk,2)) #initial value function guess
    pol_func = SharedArray{Float64}(zeros(prim.nk,2)) #initial policy function guess
    res = Results(val_func, pol_func) #initialize results struct
    prim, res #return deliverables
end

#Bellman Operator
@everywhere function Bellman(prim::Primitives,res::Results)
    @unpack val_func = res #unpack value function
    @unpack k_grid, β, δ, α, nk, z, π = prim #unpack model primitives
    v_next = SharedArray{Float64}(zeros(nk,2)) #next guess of value function to fill

    for zi = 1:2
        choice_lower = 1 #for exploiting monotonicity of policy function
        for k_index = 1:nk
            k = k_grid[k_index] #value of k
            candidate_max = -Inf #bad candidate max
            budget = z[zi]*k^α + (1-δ)*k #budget

            for kp_index in choice_lower:nk #loop over possible selections of k', exploiting monotonicity of policy function
                c = budget - k_grid[kp_index] #consumption given k' selection
                if c>0 #check for positivity
                    val = log(c) + β*(val_func[kp_index, 1]* π[zi, 1] + val_func[kp_index, 2]* π[zi, 2]) #compute value
                    if val>candidate_max #check for new max value
                        candidate_max = val #update max value
                        res.pol_func[k_index, zi] = k_grid[kp_index] #update policy function
                        choice_lower = kp_index #update lowest possible choice
                    end
                end
            end
            v_next[k_index, zi] = candidate_max #update value function
        end
    end
    v_next #return next guess of value function
end

#Value function iteration
@everywhere function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter

    while err>tol #begin iteration
        v_next = Bellman(prim, res) #spit out new vectors
        err = abs.(maximum(v_next[:,1].-res.val_func[:,1]))/abs(v_next[prim.nk, 1]) + abs.(maximum(v_next[:,2].-res.val_func[:,2]))/abs(v_next[prim.nk, 2])
        #reset error level
        res.val_func = v_next #update value function
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end

#solve the model
@everywhere function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
end

# this function will run the model in fortran
@everywhere function run_in_fortran()
    path = joinpath(pwd(), "pset1/model.f90")
    #run(`SHELL COMMANDS THAT COMPILE`)
    #run(`SHELL COMMANDS THAT RUN`)
end
