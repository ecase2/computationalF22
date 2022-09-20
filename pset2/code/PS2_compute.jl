using Parameters, Plots, Statistics #import the libraries we want
include("PS2_model.jl") #import the functions that solve our growth model

prim, res = Initialize(0.994) #initialize primitive and results structs

d = ClearMarket(prim, res)



