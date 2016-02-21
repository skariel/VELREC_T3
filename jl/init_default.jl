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

# just a utility array
c = SharedArray(eltype(rho), length(m));

to_rho!(pos,m, rho);

_s_pos = SharedArray(Float32, size(pos))
_s_vx = Array(Float32, length(vx))
_s_vy = Array(Float32, length(vx))
_s_vz = Array(Float32, length(vx))

function push_posv()
    global _s_pos
    global _s_vx
    global _s_vy
    global _s_vz
    global _s_m
    global _s_rho
    copy_into!(_s_pos, pos)
    copy_into!(_s_vx, vx)
    copy_into!(_s_vy, vy)
    copy_into!(_s_vz, vz)
    nothing
end

function pop_posv()
    copy_into!(pos, _s_pos)
    copy_into!(vx.s[:],_s_vx)
    copy_into!(vy.s[:],_s_vy)
    copy_into!(vz.s[:],_s_vz)
    nothing
end

push_posv()
