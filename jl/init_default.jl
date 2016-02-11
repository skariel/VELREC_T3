addprocs(8)

include("all.jl")

pos,vx,vy,vz,m = load_realization();
m.s[:] *= MASS_IN_BOX/sum(m);

rho = SharedArray(Complex64, (BOX_N,BOX_N,BOX_N));
c = SharedArray(eltype(rho), length(m));

to_rho!(pos,m, rho);
