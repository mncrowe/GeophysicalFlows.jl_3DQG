"""
Creates a function to run an SQG modon simulation in a finite depth
3D QG model.

Run using:

SQG_modon(dev; Nx, Ny, Nz, Lx, Ly, U, a, β, N², H, T, Ns)

Plot using:

using Plots

ψ = prob.vars.ψ
heatmap(to_CPU(ψ[:,:,Nz]'))

"""

include("3D_QG.jl")

using QGDipoles, NetCDF, Plots

function SQG_modon(dev;
                   Nx = 128,
                   Ny = 128,
                   Nz = 11,
                   Lx = 10,
                   Ly = 10,
                    U = 1.0,
                    a = 1.0,
                    β = 0.0,
                   N² = 1.0,
                    H = 1.0,
                    T = 10,
                   Ns = 10)

    if U < 0; β₀ = 0.0; else; β₀ = β; end    # values of β used in IC
    c = 0.0                                  # frame speed, (U)
    r_c, σ = Lx/2-4, 1                       # cutting parameters

    savename = "data/U_" * sstring(U) * "_a_" * sstring(a) *
               "_beta_" * sstring(β) * "_H_" * sstring(H)

    Δt = 0.5 * ((Lx/Nx)+(Ly/Ny)) / (10*abs(U))

    # define helper functions

    function mod_dist(x, y, x_c, y_c)
        return @. sqrt((mod(x + Lx/2 - x_c, Lx) - Lx/2)^2 +
            (mod(y + Ly/2 - y_c, Ly) - Ly/2)^2)
    end

    function f_cut(x, y, x_c, y_c)
        r = mod_dist(x, y, x_c, y_c)
        mask = r .< r_c
        return @. mask + exp(-(r - r_c)^2 / σ^2) * !mask
    end

    # define problem

    @info "creating problem ..."

    prob = Problem(Nz, dev;
                        nx = Nx,
                        ny = Ny,
                        Lx = Lx,
                        Ly = Ly,
                         β,
                        N²,                               
                         H,
                         U = -c,
                        dt = Δt,
                   stepper = "FilteredRK4")

    # Create modon

    @info "creating vortex IC..."

    vortex = DefSQGVortex(
        prob.grid;
        U,
        ℓ = a,
        R = [sqrt(N²)*H, Inf],
        β = β₀,
        M = 12,
        tol = 1e-8,
        K₀ = [4],
        CalcVelocity = false,
        CalcVorticity = false,
        CalcEnergy = false)

    @info "setting IC..."

    q₀ = β₀ / U * Eval_ψ_SQG(prob.grid, vortex.ψ; z = prob.params.z, U, R = [sqrt(N²)*H, Inf], β = β₀)

    q₀[:, :, 1] .= 0			# b = 0 on bottom
    q₀[:, :, Nz] = vortex.b		        # b = b_modon on top

    set_q!(prob, q₀)

    x, y = gridpoints(prob.grid)

    @info "creating save files..."

    filename = savename * ".nc"
    if isfile(filename); rm(filename); end

    nccreate(filename, "psi", "x", prob.grid.x, "y", prob.grid.y, "z", prob.params.z, "t", LinRange(0,T,Ns+1))
    nccreate(filename, "b", "x", prob.grid.x, "y", prob.grid.y, "t", LinRange(0,T,Ns+1))
    ncputatt(filename," ", Dict("H" => H, "U" => U, "a" => a, "beta" => β))

    save_field_data(prob, filename, 0)

    @info "running timestepper..."

    N_iter = T / Δt

    for i = 1:Ns

        stepforward!(prob, Int(N_iter / Ns))
        updatevars!(prob)

        @info " - iteration: " * string(Int(N_iter / Ns * i)) * ", t = " * sstring(prob.clock.t)
        if maximum(isnan.(prob.sol)); @warn "NaN detected."; end

        x_c, y_c = get_centre(prob)
        q₀ = prob.vars.q .* f_cut(x, y, x_c, y_c)

        set_q!(prob, q₀)

        save_field_data(prob, filename, Int(i))

    end

    #x_c, y_c = get_centre(prob)
    #q₀ = prob.vars.q .* f_cut(x, y, x_c, y_c)
    #heatmap(to_CPU(prob.vars.ψ[:,:,Nz]'))

end

## Functions:

to_CPU(f) = device_array(CPU())(f)

sstring(num) = string(round(num, sigdigits=2))

function save_field_data(problem, filename, i)

    grid = problem.grid
    Nx, Ny, Nz = grid.nx, grid.ny, problem.params.nz
    ψ, q = reshape(to_CPU(problem.vars.ψ), (Nx, Ny, Nz, 1)), to_CPU(problem.vars.q)

    ncwrite(ψ, filename, "psi", start = [1, 1, 1, i+1], count = [Nx, Ny, Nz, 1])
    ncwrite(q[:,:,Nz], filename, "b", start = [1, 1, i+1], count = [Nx, Ny, 1])

    return nothing

end

function get_centre(prob)

    Nz = prob.params.nz
    Lx, Ly = prob.grid.Lx, prob.grid.Ly
    x, y = prob.grid.x, prob.grid.y

    _, imax = findmax(prob.vars.q[:, :, Nz])
    _, imin = findmin(prob.vars.q[:, :, Nz])

    x_c, y_c = (x[imax[1]] + x[imin[1]]) / 2, (y[imax[2]] + y[imin[2]]) / 2

    if abs(x[imax[1]] - x[imin[1]]) > Lx / 2
        x_c = x_c < 0 ? x_c + Lx/2 : x_c - Lx/2
    end

    if abs(y[imax[2]] - y[imin[2]]) > Ly / 2
        y_c = y_c < 0 ? y_c + Ly/2 : y_c - Ly/2
    end

    return x_c, y_c
end

nothing