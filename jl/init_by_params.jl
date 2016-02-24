# REALIZATION_TYPE = "small" or "full"
# PROCS_NUM        = some integer
# MASS_REAL        = true or false
# OPT_TYPE         = "fwd", "laminar"
# SIM_TYPE         = "zel", "2lpt" or "dyn"
# ICS              = "target" or "random"
# LOGGING_FOLDER   = folder name to create / overwite

addprocs(PROCS_NUM)

include("all.jl")

try
    mkdir( "realization_$(LOGGING_FOLDER)")
catch e
end

Logging.configure(filename= "realization_$(LOGGING_FOLDER)/log.txt")

info("PROCS_NUM=",PROCS_NUM)
info("REALIZATION_TYPE=",REALIZATION_TYPE)
info("MASS_REAL=",MASS_REAL)
info("OPT_TYPE=",OPT_TYPE)
info("SIM_TYPE=",SIM_TYPE)
info("ICS=",ICS)

pos,vx,vy,vz,m = load_realization(REALIZATION_TYPE);
info("N=",length(m), " min(m)=",minimum(m)," max(m)=", maximum(m))

info("rescaling masses")

if MASS_REAL
    m.s[:] *= MASS_IN_BOX/sum(m);
else
    m.s[:] = MASS_IN_BOX/length(m);
end

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


# OPT_TYPE         = "fwd" or "laminar"
# SIM_TYPE         = "zel", "2lpt" or "dyn"
# ICS              = "target" or "random"


opos_i = SharedArray(Float32, size(pos))
if ICS == "target"
    copy_into!(opos_i, pos);
elseif ICS == "random"
    opos_i.s[:,:] = SIDE_LEN*rand(Float32, size(opos_i));
else
    throw("BAD ICS PARAMETER: "*string(ICS))
end


if     OPT_TYPE == "fwd" && SIM_TYPE == "zel"
    optimize_zeld_vs_pushed_pos!(rho, opos_i, pos, m, 0.01, 1.0, 0.21, 350);
elseif OPT_TYPE == "fwd" && SIM_TYPE == "2lpt"
    optimize_2lpt_vs_pushed_pos!(rho, opos_i, pos, m, 0.01, 1.0, 0.21, 350);
elseif OPT_TYPE == "fwd" && SIM_TYPE == "dyn"
    optimize_dyn_vs_pushed_pos!(vx,vy,vz, rho, opos_i, pos, m, 0.01, 1.0, 0.21, 350);
elseif OPT_TYPE == "laminar" && SIM_TYPE == "zel"
    back_optimize_zeld_vs_pushed_pos!(rho, opos_i, pos, m, 0.01, 1.0, 50, 0.21, 350);
elseif OPT_TYPE == "laminar" && SIM_TYPE == "2lpt"
    back_optimize_2lpt_vs_pushed_pos!(rho, opos_i, pos, m, 0.01, 1.0, 50, 0.21, 350);
elseif OPT_TYPE == "laminar" && SIM_TYPE == "dyn"
    back_optimize_dyn_vs_pushed_pos!(vx,vy,vz, rho, opos_i, pos, m, 0.01, 1.0, 50, 0.21, 350);
else
    throw("BAD OPT_TYPE AND/OR SIM_TYPE PARAMETERS: "*string(OPT_TYPE)+", "*string(SIM_TYPE))
end

save_realization(opos_i,vx,vy,vz,m,LOGGING_FOLDER)
