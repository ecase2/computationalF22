#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 7
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#
using Parameters, Distributions, Random, Optim, PyPlot, Plots, LinearAlgebra

# Structure of model primitives (unmutable objects)
@with_kw struct Params
    ### true process: xt = ρ0*x(t-1) + ε, where ε ∼ N(0, σ0^2).
    ### parameters to be estimated are b = (ρ0, σ0^2)

    ρ0::Float64 = 0.5 # autocorrelation of true AR(1) process
    σ0::Float64 = 1.0 # standard deviation of true AR(1) process
    y0::Float64 = 0.0 # the first element of true AR(1) process
    T::Int64 = 200 # the length of x series
    l::Int64 = 2 # number of parameters to be estimated
    n::Int64 = 3 # number of moments to be used for estimation

    H::Int64 = 10

    type::Int64 = 1 ### type stands for exercises:
    ### type 1 = ex.4 - use mean and variance
    ### type 2 = ex.5 - use variance and covariance
    ### type 3 = ex.6 - use mean, variance and covariance

    iT::Int64 = 4

end

mutable struct Results
    # for true data
    ε::Array{Float64}
    y::Array{Float64}

    # for simulated data
    ε_grid::Array{Float64, 2}
    y_grid::Array{Float64, 2}

    W::Array{Float64, 2}

    # store estimates 1 and 2
    # b11 = ρ estimate with W = I (=ρ1 in the handout)
    # b21 = σ estimate with W = I (=σ1 in the handout)
    # b12 = ρ estimate with W = to get efficient estimate (=ρ2 in the handout)
    # b22 = σ estimate with W = to get efficient estimate (=σ2 in the handout)

    b11::Float64
    b21::Float64

    b12::Float64
    b22::Float64

end

# Function to initialize parameters and results structures
function Initialize(; type::Int64 = 1)

    if type == 1 || type == 2
        par = Params(n = 2, type = type)
    elseif type == 3
        par = Params(n = 3, type = type)
    end


    ε = zeros(par.T)
    y = zeros(par.T)

    ε_grid = zeros(par.T, par.H)
    y_grid = zeros(par.T, par.H)

    W = zeros(par.n, par.n)

    b11, b21, b12, b22 = 0, 0, 0, 0

    res = Results(ε, y, ε_grid, y_grid, W, b11, b21, b12, b22)

    return par, res
end


function true_data(par::Params, res::Results)
    @unpack ρ0, σ0, y0, T = par
    @unpack ε, y = res

    distr = Normal(0, σ0)
    Random.seed!(123)
    ε = rand(distr, T)

    y[1] = y0
    for i = 2:T
        y[i] = ρ0*y[i-1] + ε[i]
    end
    return y

end

function errors_sim(par::Params, res::Results)
    @unpack T, H = par
    @unpack ε_grid = res

    distr = Normal(0, 1)
    Random.seed!(123)
    ε_grid = reshape(rand(distr, T*H), T, H)

    return ε_grid

end


function moments(par::Params, res::Results)
    @unpack T, H, n = par
    @unpack y, y_grid = res

    M_data = zeros(n)
    M_simulations = zero(n)

    mean_data = 1/T*sum(y)
    mean_sim  = 1/(T*H)*sum(y_grid)

    var_data = 1/T*sum((y .- mean_data).^2)
    var_sim = 1/(T*H)*sum((y_grid .- mean_data).^2)

    diff_data = y .- mean_data
    diff_sim = y_grid .- mean_sim

    acov_data, acov_sim = 0, 0

    for i = 2:T
        acov_data = acov_data + diff_data[i-1]*diff_data[i]
    end
    for i = 2:T
        for j = 1:H
            acov_sim = acov_sim + diff_sim[i-1, j]*diff_sim[i, j]
        end
    end
    acov_data = 1/T*acov_data
    acov_sim = 1/(T*H)*acov_sim

    return mean_data, mean_sim, var_data, var_sim, acov_data, acov_sim

end

function J(x::Vector)
    @unpack T, H, type, y0 = par
    @unpack ε_grid, y_grid, ε, y, W = res

    # simulate data
    # x[1] = ρ, x[2] = σ

    # not sure about the first element in the simulated matrix
    y_grid[1, :] .= 0

    for i = 2:T, j = 1:H
        y_grid[i, j] = x[1]*y_grid[i-1, j] + x[2]*ε_grid[i, j]
    end
#    println("y_grid = $y_grid\n")
#    println("ε_grid = $ε_grid")

    mean_data, mean_sim, var_data, var_sim, acov_data, acov_sim = moments(par::Params, res::Results)

    if type == 1

        M_data, M_sim = zeros(2), zeros(2)
        M_data[1], M_data[2] = mean_data, var_data
        M_sim[1], M_sim[2] = mean_sim, var_sim

    elseif type == 2

        M_data, M_sim = zeros(2), zeros(2)
        M_data[1], M_data[2] = var_data, acov_data
        M_sim[1], M_sim[2] = var_sim, acov_sim

    elseif type == 3

        M_data, M_sim = zeros(3), zeros(3)
        M_data[1], M_data[2], M_data[3] = mean_data, var_data, acov_data
        M_sim[1], M_sim[2], M_sim[3] = mean_sim, var_sim, acov_sim

    end

    return transpose(M_data - M_sim)*W*(M_data - M_sim)
end

function dif_for_gamma(par::Params, res::Results)
    @unpack T, H, type, n, iT = par
    @unpack b11, b21, y_grid, ε_grid = res

    y_grid[1, :] .= 0

    for i = 2:T, j = 1:H
        y_grid[i, j] = b11*y_grid[i-1, j] + b21*ε_grid[i, j]
    end

    mean_data, mean_sim, var_data, var_sim, acov_data, acov_sim = moments(par, res)
    if type == 1
        M_TH = zeros(n)
        M_TH[1], M_TH[2] = mean_sim, var_sim
        m = zeros(n, T, H)
        dif = zeros(n, T, H)
        m[1, :, :] = y_grid
        m[2, :, :] = (y_grid .- M_TH[1]).^2
        dif[1, :, :] = m[1, :, :] .- M_TH[1]
        dif[2, :, :] = m[2, :, :] .- M_TH[2]
    end

    if type == 2
        M_TH = zeros(2)
        M_TH[1], M_TH[2] = var_sim, acov_sim
        mean_sim = sum(y_grid)/(T*H)
        m = zeros(n, T, H)

        dif = zeros(n, T, H)

        m[1, :, :] = (y_grid .- M_TH[1]).^2
        d = y_grid .- mean_sim
        for i = 2:T
            for j = 1:H
                m[2, i, j] = d[i-1, j]*d[i, j]
            end
        end
        dif[1, :, :] = m[1, :, :] .- M_TH[1]
        dif[2, :, :] = m[2, :, :] .- M_TH[2]

    end

    if type == 3
        M_TH = zeros(3)
        M_TH[1], M_TH[2], M_TH[3] = mean_sim, var_sim, acov_sim

        m = zeros(n, T, H)
        dif = zeros(n, T, H)
        m[1, :, :] = y_grid
        m[2, :, :] = (y_grid .- M_TH[1]).^2

        d = y_grid .- M_TH[1]
        for i = 2:T
            for j = 1:H
                m[3, i, j] = d[i-1, j]*d[i, j]
            end
        end
        dif[1, :, :] = m[1, :, :] .- M_TH[1]
        dif[2, :, :] = m[2, :, :] .- M_TH[2]
        dif[3, :, :] = m[3, :, :] .- M_TH[3]

    end
    return dif

end

function gamma(par::Params, res::Results; dif::Array{Float64, 3})
    @unpack T, H, iT, n = par

    Γ0 = zeros(n, n)
    for t = 1:T, h = 1:H
        v = dif[:, t, h]
        Γ0 = Γ0 + v*transpose(v)
    end
    Γ0 = Γ0/(T*H)


    Γj = zeros(n, n, iT)
    for j = 1:iT, h = 1:H, t = j + 1
        v1 = dif[:, t, h]
        v2 = dif[:, 1, h]
        Γj[:, :, j] += v1*transpose(v2)
    end
    Γj = Γj/(T*H)

    return Γ0, Γj
end

function S(par::Params, res::Results; Γ0::Matrix{Float64}, Γj::Array{Float64, 3})
    @unpack iT, H, n = par

    s = Γ0

    for j = 1:iT
        s = s + (1- j/(iT+1))*(Γj[:, :, j] + transpose(Γj[:, :, j]))
    end

    S_hat = (1+1/H)*s

    return S_hat
end



### PART FOR standard errors - use the N-W estimates
### Computation of derivatives

function model_data(par::Params, res::Results; b1::Float64, b2::Float64)
    @unpack T, H = par
    @unpack y_grid, ε_grid = res

    y_grid[1, :] .= 0
#    println("b1 = $b1, b2 = $b2")
    for i = 2:T, j = 1:H
        y_grid[i, j] = b1*y_grid[i-1, j] + b2*ε_grid[i, j]
    end

    return y_grid

end

function derivatives(par::Params, res::Results; s::Float64 = 0.01)
    @unpack type = par
    @unpack ε_grid, y_grid, b12, b22 = res

    res.y_grid = model_data(par,res; b1 = b12, b2 = b22)
#    println("y_grid = ", res.y_grid[3, 1:10], "\n")
    m_d, m_sim, var_d, var_sim, acov_d, acov_sim = moments(par::Params, res::Results)
#    println("Moments 1: m = $m_sim, var = $var_sim, acov = $m_sim \n")

    res.y_grid = model_data(par,res; b1 = b12 - s, b2 = b22)
#    println("y_grid = ", res.y_grid[3, 1:10], "\n")
    m_d_ρ, m_sim_ρ, var_d_ρ, var_sim_ρ, acov_d_ρ, acov_sim_ρ = moments(par::Params, res::Results)
#    println("Moments 2: m  = $m_sim_ρ, var = $var_d_ρ, acov = $m_sim_ρ \n")
    der_ρ = [-(m_sim - m_sim_ρ)/s, -(var_sim - var_sim_ρ)/s, -(acov_sim - acov_sim_ρ)/s]
#    println("1 der $der_ρ")

    res.y_grid = model_data(par,res; b1 = b12, b2 = b22 - s)
#    println("y_grid = ", res.y_grid[3, 1:10], "\n")
    m_d_σ, m_sim_σ, var_d_σ, var_sim_σ, acov_d_σ, acov_sim_σ = moments(par::Params, res::Results)
#    println("Moments 3: m = $m_sim_σ, var = $var_sim_σ, acov =  $m_sim_σ \n")
    der_σ = [-(m_sim - m_sim_σ)/s, -(var_sim - var_sim_σ)/s, -(acov_sim - acov_sim_σ)/s]
#    println("1 der $der_σ")

    if type == 1
        der = [der_ρ[1] der_σ[1];
                der_ρ[2] der_σ[2]]
    elseif type == 2
        der = [der_ρ[2] der_σ[2];
                der_ρ[3] der_σ[3]]
    elseif type == 3
        der = [der_ρ[1] der_σ[1];
                der_ρ[2] der_σ[2];
                der_ρ[3] der_σ[3]]
    end

    return der

end

function st_errors(par::Params; ∇::Matrix{Float64}, S_hat::Matrix{Float64})
    @unpack T, type = par

    var_cov = 1/T*inv(transpose(∇)*inv(S_hat)*∇)

    return sqrt(var_cov[1, 1]), sqrt(var_cov[2, 2])
end


function J_test(par::Params; J::Float64 = sol2.minimum)
    @unpack T, H = par

    return T*H/(1+H)*J
end

######## FOR PLOTS
function J_func(par::Params, res::Results; ρ::Float64, σ::Float64)
    @unpack T, H, type = par
    @unpack ε_grid, y_grid, ε, y, W = res

    ρ = ρ
    σ = σ

    y_grid[1, :] .= 0

    for i = 2:T, j = 1:H
        y_grid[i, j] = ρ*y_grid[i-1, j] + σ*ε_grid[i, j]
    end

    mean_data, mean_sim, var_data, var_sim, acov_data, acov_sim = moments(par::Params, res::Results)

    if type == 1

        M_data, M_sim = zeros(2), zeros(2)
        M_data[1], M_data[2] = mean_data, var_data
        M_sim[1], M_sim[2] = mean_sim, var_sim

    elseif type == 2

        M_data, M_sim = zeros(2), zeros(2)
        M_data[1], M_data[2] = var_data, acov_data
        M_sim[1], M_sim[2] = var_sim, acov_sim

    elseif type == 3

        M_data, M_sim = zeros(3), zeros(3)
        M_data[1], M_data[2], M_data[3] = mean_data, var_data, acov_data
        M_sim[1], M_sim[2], M_sim[3] = mean_sim, var_sim, acov_sim

    end

    return transpose(M_data - M_sim)*W*(M_data - M_sim)
end

function plot_J(par::Params, res::Results; ρl::Float64 = 0.35, ρstep::Float64 = 0.05, ρh::Float64 = 0.65, σl::Float64 = 0.8, σstep::Float64 = 0.05,
    σh::Float64 = 1.2)

    ρ_grid = collect(ρl:ρstep:ρh)
    σ_grid = collect(σl:σstep:σh)

    nρ = length(ρ_grid)
    nσ = length(σ_grid)

    z_ax = zeros(nρ, nσ)
    x_ax = zeros(nρ, nσ)
    y_ax = zeros(nρ, nσ)

    for i = 1:nρ, j = 1:nσ
        z_ax[i, j] = J_func(par, res; ρ = ρ_grid[i], σ = σ_grid[j])
        x_ax[i, j] = ρ_grid[i]
        y_ax[i, j] = σ_grid[j]
    end

    x_ax = vec(x_ax)
    y_ax = vec(y_ax)
    z_ax = vec(z_ax)

    return x_ax, y_ax, z_ax
end

function solve_model(; type::Int64 = 1)
    par, res = Initialize(; type)
    res.y = true_data(par, res) # it is not dependent on the type (case) and will be the same for all types.

    # Task 3. Sumulate errors for simulations.
    res.ε_grid = errors_sim(par, res) # will be the same for all types

    # Task 4.
    # (a) Plot and consistent estimate b1 with W = I
    res.W = Diagonal(ones(par.n))

    # Plot
    x_ax, y_ax, z_ax = plot_J(par, res; ρl = 0.35, ρstep = 0.05, ρh = 0.65, σl = 0.8, σstep = 0.05, σh = 1.2)
    Plots.plot(x_ax, y_ax, z_ax, st=:surface, camera=(-15, 15))
    # add saving option

    # Estimates
    sol = optimize(J, [0.5, 1.0])
    res.b11, res.b21 = sol.minimizer[1], sol.minimizer[2]

    # (b) Find W* and efficient estimate b2 with W = W*
    dif = dif_for_gamma(par, res)
    Γ0, Γj = gamma(par, res; dif = dif)
    S_hat = S(par, res; Γ0 = Γ0, Γj = Γj)
    res.W = inv(S_hat)

    # Estimates
    sol2 = optimize(J, [0.5, 1.0])
    res.b12, res.b22 = sol2.minimizer[1], sol2.minimizer[2]

    # (c) Derivatives and standard errors
    ∇ = derivatives(par, res)
    errors = st_errors(par; ∇ = ∇, S_hat = S_hat)

    # (d) J-test WRONG
    j_st = J_test(par)

    return res.b11, res.b12, res.W, ∇, errors, j_st

end
