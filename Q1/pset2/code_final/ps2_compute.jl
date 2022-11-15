using Parameters, Plots, Statistics #import the libraries we want
include("PS2_model.jl") #import the functions that solve our growth model

prim, res = Initialize(0.994) #initialize primitive and results structs

### Solve the model 
### Return net wealth
d = ClearMarket(prim, res)

### Compute the Gini coefficient
sum_wealth, sum_people = Lorenz(prim, res)
Gini_coef = Gini(prim, res)

# make graphs and lorenz curve
include("ps2_graphs.jl")

