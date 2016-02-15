addprocs(8)

include("all.jl")

pos,vx,vy,vz,m = load_realization("small");
info("N=",length(m), " min(m)=",minimum(m)," max(m)=", maximum(m))

#pos,vx,vy,vz,m = filter_mass_above(pos,vx,vy,vz,m,310.0)
# info("N=",length(m), " min(m)=",minimum(m)," max(m)=", maximum(m))

info("rescaling masses")
m.s[:] *= MASS_IN_BOX/sum(m);
#m.s[:] = MASS_IN_BOX/length(m);
info("N=",length(m), " min(m)=",minimum(m)," max(m)=", maximum(m))

rho = SharedArray(Complex64, (BOX_N,BOX_N,BOX_N));
c = SharedArray(eltype(rho), length(m));

to_rho!(pos,m, rho);
