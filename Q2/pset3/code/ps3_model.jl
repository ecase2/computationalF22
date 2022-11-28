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
