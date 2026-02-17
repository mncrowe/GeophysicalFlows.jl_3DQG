include("SQG_modon.jl")

# set parameters

dev = GPU()
Nx, Ny, Nz = 2048, 2048, 13
#Nx, Ny, Nz = 512, 512, 9
Lx, Ly = 40.96, 40.96

U, a, β = -1.0, 1.0, 0.4
N², H = 1.0, 1.0	

T = 50
Ns = 100

SQG_modon(dev; Nx, Ny, Nz, Lx, Ly, U, a, β, N², H, T, Ns)