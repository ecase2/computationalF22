using Parameters, Plots #import the libraries we want


#keyword-enabled structure to hold model primitives
@with_kw struct Primitives
    β::Float64 = 0.9932 # discount rate
    α::Float64 = 1.5 # relative risk-aversion
    a_min::Float64 = -2 # assets lower bound
    a_max::Float64 = 5.0 # assets upper bound
    na::Int64 = 400 #number of assets grid points
    a_grid::Array{Float64,1} = collect(range(a_min, length = na, stop = a_max)) # assets grid
    ns::Int64 = 2 # number of possible statuses (employed or unemployed)
    s_grid::Array{Float64, 1} = [1.0, 0.5] # earnings when employed and unemployed
    trans_matrix = [0.97 0.03; 0.50 0.50]
end

#structure that holds model results
# With stocahstic z the objects become two-dimensional
mutable struct Results
    val_func::Array{Float64, 2} # value function
    pol_func::Array{Float64, 2} # policy function
    distr::Array{Float64, 2} # distribution
    q::Float64 # price of a non state-contingent bond
    λ::Array{Float64, 2}
end

#function for initializing model primitives and results
function Initialize(qval)
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.na, prim.ns) #initial value function guess
    pol_func = zeros(prim.na, prim.ns) #initial policy function guess
    distr = zeros(prim.na, prim.ns)
    q = qval 
    lgrid = zeros(prim.na, prim.ns)
    res = Results(val_func, pol_func, distr, q, lgrid) #initialize results struct
    prim, res #return deliverables
end

#Bellman Operator
function Bellman(prim::Primitives, res::Results)
    @unpack val_func, q, pol_func = res #unpack value function
    @unpack a_grid, β, α, ns, s_grid, trans_matrix, na  = prim #unpack model primitives
    v_next = zeros(na, ns) #next guess of value function to fill

    for s_index = 1:ns
        s = s_grid[s_index] #value of s - current status
        choice_lower = 1 #for exploiting monotonicity of policy function

        for a_index = 1:na
            candidate_max = -Inf

            a = a_grid[a_index] #value of a - current assets level

            budget = s + a #budget
            for ap_index in choice_lower:na #loop over possible selections of k', exploiting monotonicity of policy function
                c = budget - q*a_grid[ap_index] #consumption given k' selection
                if c>=0 #check for positivity
                    util = (c^(1-α)-1)/(1-α)
                    val = util + β*transpose(trans_matrix[s_index, :])*val_func[ap_index, :]   #compute value
                    if val>candidate_max #check for new max value
                        candidate_max = val #update max value
                        res.pol_func[a_index, s_index] = a_grid[ap_index] #update policy function
                        choice_lower = ap_index #update lowest possible choice
                    end
                end 
            end
            v_next[a_index, s_index] = candidate_max #update value function
        end
    end
    v_next, res.pol_func #return next guess of value function
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    @unpack pol_func, val_func = res
    n = 0 #counter
    
    while err>tol #begin iteration
        v_next, res.pol_func = Bellman(prim, res) #spit out new vectors
        err = abs.(maximum(v_next.-res.val_func))/abs(v_next[prim.na, 1]) #reset error level
        res.val_func = v_next #update value function
        n+=1
        #println("pol_func[1:10, 1] = ", res.pol_func[1:10, 1])
    end
#    println("Value function converged in ", n, " iterations.")
    return pol_func 
end

#solve the model
#function Solve_model(prim::Primitives, res::Results)
#    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
    #return res.pol_func
#end

### Indicator matrix - to know whether the agent with (a,s) choose a'.
function Indicator(prim::Primitives, res::Results)
    @unpack val_func, q, pol_func = res #unpack value function
    @unpack a_grid, β, α, ns, s_grid, trans_matrix, na = prim #unpack model primitives
    ### create the indicator variable
    I = zeros(na, na, ns)
    for s_index in 1:ns
        for a_index in 1:na
            for ap_index in 1:na
                ap = a_grid[ap_index]
                if res.pol_func[a_index, s_index] == ap
                    I[a_index, ap_index, s_index] = 1
                end
            end
        end
    end
    return I
end

### Function to find stationary distribution
function Distr(prim::Primitives, res::Results; iter = 1000, tol = 0.000001)
    @unpack val_func, q, pol_func, distr = res #unpack value function
    @unpack a_grid, β, α, ns, s_grid, trans_matrix, na = prim #unpack model primitives
    I = Indicator(prim, res)
    temp_distr = ones(na, ns)/(na*ns)
    ### Update intial guess and difference according to Dean's notes
    diff = 100
    n = 1
    ### create the indicator variable
    while (diff > tol && n < iter)
        distr_1 = zeros(na, ns)
        for sp_index in 1:ns
            for ap_index in 1:na
#                distr_1[ap_index, sp_index] = distr_1[ap_index, sp_index] + sum(trans_matrix[1, sp_index].*(transpose(temp_distr[:, 1]).*transpose(I[:, ap_index, sp_index])) +  trans_matrix[2, sp_index].*(transpose(temp_distr[:, 2]).*transpose(I[:, ap_index, sp_index])))
                distr_1[ap_index, sp_index] = distr_1[ap_index, sp_index] + sum(trans_matrix[1, sp_index].*(transpose(temp_distr[:, 1]).*transpose(I[:, ap_index, sp_index])) +
                 trans_matrix[2, sp_index].*(transpose(temp_distr[:, 2]).*transpose(I[:, ap_index, sp_index])))
### NOT SURE that 1st is correct
            end
        end
    #    diff = sqrt(sum((distr_1 .- temp_distr).^2))
        diff = maximum(abs.(distr_1 .- temp_distr))
        n = n + 1
    #    println("Iteration ", n-1, " Diff = ", diff);
        temp_distr = distr_1
    end
    res.distr = temp_distr
    #return temp_distr
end

function ExcessDemand(prim::Primitives, res::Results)
#    @unpack val_func, q, pol_func, distr = res #unpack value function
#    @unpack a_grid, β, α, ns, s_grid, trans_matrix, na = prim #unpack model primitives

    ED = transpose(res.distr[:, 1])*res.pol_func[:, 1] + transpose(res.distr[:, 2])*res.pol_func[:, 2]
    return ED
end

function ClearMarket(prim::Primitives, res::Results; iter = 1000, tol = 0.005)
    @unpack val_func, q, pol_func, distr = res #unpack value function
    @unpack a_grid, β, α, ns, s_grid, trans_matrix, na = prim #unpack model primitives
    V_iterate(prim, res)
    Distr(prim, res)
    diff = ExcessDemand(prim, res)
    n = 1
    while (abs(diff) > tol && n < iter)
        adj_step = 0.00005
        if abs(diff) <0.06
            adj_step = 0.000001
        end 
        res.q = res.q + adj_step
        Initialize(res.q)
        V_iterate(prim, res)
        Distr(prim, res)
        diff = ExcessDemand(prim, res)
        println("Iteration ", n-1, " Diff = ", diff, " under price = ", res.q);
        n = n + 1
    end
    return diff

end

# Calculate first-best welfare
function calcWelfareFB(prim::Primitives)
    @unpack β, α, trans_matrix = prim #unpack model primitives
    π = trans_matrix
    # Calculate the average fraction of agents who are employed
    U_bar = π[1,2]/(1 - π[2,2] + π[1,2])
    E_bar = 1 - U_bar 

    # Calculate consumption
    c_bar = E_bar + 0.5*(1-E_bar)

    # Calculate FB welfare
    W_FB = ((c_bar^(1-α) - 1)/(1-α))/(1-β)

    # Return FB value
    W_FB
end

# Finds consumption equivalents
function calcCE(prim::Primitives, res::Results)
    @unpack val_func, λ = res #unpack value function and consumption equivalence
    @unpack β, α, na, ns = prim #unpack model primitives

    # Calculate parameter value
    θ = 1/((1-α)*(1-β))
    W_FB = calcWelfareFB(prim)

    # Calculate consumption equivalents
    for s_index = 1:ns
        for a_index = 1:na
            λ[a_index,s_index] = ((W_FB + θ)/(val_func[a_index,s_index] + θ))^(1/(1-α)) - 1
        end
    end

    # Return λ
    λ
end

# Calculate welfare from incomplete markets
function calcWelfareInc(prim::Primitives, res::Results)
    @unpack val_func, distr = res
    @unpack na, ns = prim

    # Calculate W_Inc
    W_Inc = 0

    for s_index = 1:ns
        for a_index = 1:na
            W_Inc = W_Inc .+ distr[a_index,s_index].*val_func[a_index,s_index]
        end
    end

    # Return W_Inc
    W_Inc
end

# Calculate welfare gain
function calcWelfareGain(prim::Primitives, res::Results)
    @unpack distr, λ = res
    @unpack na, ns = prim

     # Calculate WG
     WG = 0
     
     for s_index = 1:ns
         for a_index = 1:na
             WG = WG .+ λ[a_index,s_index].*distr[a_index,s_index]
         end
     end
 
     # Return WG
     WG
end


function voteforcompmarkets(res)
    @unpack distr, λ = res

    l = zeros(size(λ))
    summ = 0.0 
    for i = 1:size(l)[1]
        for j = 1:2
            if λ[i,j] > 0.0
                l[i,j] = distr[i,j]
                summ += l[i,j]
            end
        end
    end
    summ

end







#######################################################################################
### Function to find stationary distribution when there is only one state of the world
function DistrSimple(prim::Primitives, res::Results; iter = 1000, tol = 0.000001)
    @unpack val_func, q, pol_func = res #unpack value function
    @unpack a_grid, β, α, ns, s_grid, trans_matrix, na = prim #unpack model primitives

#    temp_distr = ones(na, ns)/(na*ns)
#    distr_1 = zeros(na, ns)

    temp_distr = ones(na)/na
    diff = 100
    n = 1
    ### create the indicator variable
    while (diff > tol && n < iter)
            distr_1 = zeros(na)
    #    for sp_index in 1:ns
            for a_index in 1:na
                ap = pol_func[a_index, 1]
                ap_index = findall(x -> x ==  ap, a_grid)[1]

            #    distr_1[ap_index, sp_index] = distr_1[ap_index, sp_index] + sum(trans_matrix[1, sp_index].*(transpose(temp_distr[:, 1]).*transpose(I[:, ap_index, sp_index])) +  trans_matrix[2, sp_index].*(transpose(temp_distr[:, 2]).*transpose(I[:, ap_index, sp_index])))
                distr_1[ap_index] = distr_1[ap_index] + temp_distr[a_index]
            end

        diff = sqrt(sum((distr_1 .- temp_distr).^2))
        n = n + 1
        println("Iteration ", n-1, " Diff = ", diff);
        temp_distr = distr_1
    end
    return temp_distr
end
