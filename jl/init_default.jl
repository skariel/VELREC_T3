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

_s_pos = nothing
_s_vx = nothing
_s_vy = nothing
_s_vz = nothing
_s_m = nothing
_s_rho = nothing

function push_realization()
    global _s_pos
    global _s_vx
    global _s_vy
    global _s_vz
    global _s_m
    global _s_rho
    _s_pos = pos.s[:,:]
    _s_vx = vx.s[:]
    _s_vy = vy.s[:]
    _s_vz = vz.s[:]
    _s_m = m.s[:]
    _s_rho = rho.s[:,:,:]
    nothing
end

function pop_realization()
    pos.s[:,:] = _s_pos
    vx.s[:] = _s_vx
    vy.s[:] = _s_vy
    vz.s[:] = _s_vz
    m.s[:] = _s_m
    rho.s[:,:,:] = _s_rho
    nothing
end

push_realization()
