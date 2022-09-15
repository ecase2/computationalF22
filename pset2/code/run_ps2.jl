using Parameters
include("model_ps2.jl")

θ, res = initialize()
@time solve_model(res)