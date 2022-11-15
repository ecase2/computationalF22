#=
    PROJECT:     COMPUTATIONAL FALL 2022 PSET 7
    AUTHORS:     Hanna Han, Anna Lukianova
    CONTENTS:    This file contains main functions.
=#

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

    # true data moments
    m_data::Float64
    var_data::Float64
    acor_data::Float64

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

    m_data, var_data, acor_data = 0, 0, 0

    res = Results(ε, y, ε_grid, y_grid, W, b11, b21, b12, b22, m_data, var_data, acor_data)

    return par, res
end


function true_data(par::Params, res::Results)
    @unpack ρ0, σ0, y0, T = par

    distr = Normal(0, σ0)
    Random.seed!(123)
    res.ε = rand(distr, T)

    res.y[1] = y0
    for i = 2:T
        res.y[i] = ρ0*res.y[i-1] + res.ε[i]
    end

    res.m_data = mean(res.y)
    res.var_data = var(res.y)
    a = autocor(res.y, [1]; demean = true)
    res.acor_data = a[1]

    return res.y, res.m_data, res.var_data, res.acor_data

end

function errors_sim(par::Params, res::Results)
    @unpack T, H = par

    distr = Normal(0, 1)
    Random.seed!(123)
    res.ε_grid = reshape(rand(distr, T*H), T, H)

    return res.ε_grid

end

function sim_data(par::Params, res::Results; b1::Float64, b2::Float64)
    @unpack T, H = par

    res.y_grid[1, :] .= 0 .+ res.ε_grid[1, :]

    for i = 2:T, j = 1:H
        res.y_grid[i, j] = b1*res.y_grid[i-1, j] + b2*res.ε_grid[i, j]
    end


    m_sim = mean(res.y_grid)
    var_sim = var(res.y_grid)
    a = autocor(res.y_grid, [1]; demean = true)
    acor_sim = a[1]

    return res.y_grid, m_sim, var_sim, acor_sim

end

function J(x::Vector)
    @unpack T, H, type, y0 = par
    @unpack ε_grid, y_grid, ε, y, W, m_data, var_data, acor_data  = res

    # simulate data and find moments of the simulated data

    res.y_grid, m_sim, var_sim, acor_sim = sim_data(par::Params, res::Results; b1 = x[1], b2 = x[2])

    if type == 1

        M_data = [m_data, var_data]
        M_sim = [m_sim, var_sim]

    elseif type == 2

        M_data = [var_data, acor_data]
        M_sim = [var_sim, acor_sim]

    elseif type == 3

        M_data = [m_data, var_data, acor_data]
        M_sim = [m_sim, var_sim, acor_sim]

    end
    J_val = transpose(M_data - M_sim)*W*(M_data - M_sim)
    return J_val
end

function dif_for_gamma(par::Params, res::Results)
    @unpack T, H, type, n, iT = par
    @unpack b11, b21, y_grid, ε_grid = res

    res.y_grid, m_sim, var_sim, acor_sim = sim_data(par::Params, res::Results; b1 = b11, b2 = b21)

    if type == 1
        M_TH = [m_sim, var_sim]
        m = zeros(n, T, H)
        dif = zeros(n, T, H)
        m[1, :, :] = y_grid
        m[2, :, :] = (y_grid .- m_sim).^2
        dif[1, :, :] = m[1, :, :] .- M_TH[1]
        dif[2, :, :] = m[2, :, :] .- M_TH[2]
    end

    if type == 2
        M_TH = [var_sim, acor_sim]
        m = zeros(n, T, H)
        dif = zeros(n, T, H)
        d = y_grid .- m_sim
        m[1, :, :] = d.^2
        for i = 2:T
            for j = 1:H
                m[2, i, j] = d[i-1, j]*d[i, j]
            end
        end
        dif[1, :, :] = m[1, :, :] .- M_TH[1]
        dif[2, :, :] = m[2, :, :] .- M_TH[2]

    end

    if type == 3
        M_TH = [m_sim, var_sim, acor_sim]

        m = zeros(n, T, H)
        dif = zeros(n, T, H)
        d = y_grid .- m_sim
        m[1, :, :] = y_grid
        m[2, :, :] = d.^2

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

function derivatives(par::Params, res::Results; s::Float64 = 0.01)
    @unpack type = par
    @unpack ε_grid, y_grid, b12, b22 = res

    res.y_grid, m_sim, var_sim, acor_sim = sim_data(par,res; b1 = b12, b2 = b22)

    res.y_grid, m_sim_ρ, var_sim_ρ, acor_sim_ρ = sim_data(par,res; b1 = b12 - s, b2 = b22)
    der_ρ = [-(m_sim - m_sim_ρ)/s, -(var_sim - var_sim_ρ)/s, -(acor_sim - acor_sim_ρ)/s]

    res.y_grid, m_sim_σ, var_sim_σ, acor_sim_σ = sim_data(par,res; b1 = b12, b2 = b22 - s)
    der_σ = [-(m_sim - m_sim_σ)/s, -(var_sim - var_sim_σ)/s, -(acor_sim - acor_sim_σ)/s]

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

    return [sqrt(var_cov[1, 1]), sqrt(var_cov[2, 2])]
end


function J_test(par::Params; J::Float64 = sol2.minimum)
    @unpack T, H = par

    return T*H/(1+H)*J
end

######## FOR PLOTS
function J_func(par::Params, res::Results; ρ::Float64, σ::Float64)
    @unpack T, H, type = par
    @unpack ε_grid, y_grid, ε, y, W, m_data, var_data, acor_data = res

    ρ = ρ
    σ = σ

    res.y_grid, m_sim, var_sim, acor_sim = sim_data(par::Params, res::Results; b1 = ρ, b2 = σ)

    if type == 1

        M_data = [m_data, var_data]
        M_sim = [m_sim, var_sim]

    elseif type == 2

        M_data = [var_data, acor_data]
        M_sim = [var_sim, acor_sim]

    elseif type == 3

        M_data = [m_data, var_data, acor_data]
        M_sim = [m_sim, var_sim, acor_sim]

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


function bootstrap(par::Params, res::Results; N::Int64 = 1000)
    # N is number of simulations
        @unpack T, H, σ0, y0, ρ0 = par
        @unpack b11, b21 = res

        b1 = zeros(N, 2)
        b2 = zeros(N, 2)


        for i = 1:N


            # True data
            distr = Normal(0, σ0)
            res.ε = rand(distr, T)

            res.y[1] = y0
            for i = 2:T
                res.y[i] = ρ0*res.y[i-1] + res.ε[i]
            end

            res.m_data = mean(res.y)
            res.var_data = var(res.y)
            a = autocor(res.y, [1]; demean = true)
            res.acor_data = a[1]

            # Simulated data

            distr1 = Normal(0, 1)
            res.ε_grid = reshape(rand(distr1, T*H), T, H)

            res.W = Diagonal(ones(par.n))
            sol = optimize(J, [0.5, 1.0])
            res.b11, res.b21 = sol.minimizer[1], sol.minimizer[2]

            # (b) Find W* and efficient estimate b2 with W = W*
            dif = dif_for_gamma(par, res)
            Γ0, Γj = gamma(par, res; dif = dif)
            S_hat = S(par, res; Γ0 = Γ0, Γj = Γj)
            res.W = inv(S_hat)

            # Estimates
            sol = optimize(J, [0.5, 1.0])
            b12, b22 = sol.minimizer[1], sol.minimizer[2]

            b1[i, :] = [res.b11, res.b21]
            b2[i, :] = [b12, b22]
        end

    return b1, b2
end
