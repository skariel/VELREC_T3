# RELNAME
# PROCS_NUM
# A_FROM
# A_TO
# LOGGING_FOLDER


ICORE = true

addprocs(PROCS_NUM)

include("../all.jl")

try
    mkdir( "realization_$(LOGGING_FOLDER)")
catch e
end

Logging.configure(filename= "realization_$(LOGGING_FOLDER)/log.txt")

info("PROCS_NUM=",PROCS_NUM)
info("RELNAME=",RELNAME)
info("A_FROM=",A_FROM)
info("A_TO=",A_TO)

pos,vx,vy,vz,m = load_realization(RELNAME);
info("N=",length(m), " min(m)=",minimum(m)," max(m)=", maximum(m))

rho = SharedArray(Complex64, (BOX_N,BOX_N,BOX_N));

# just a utility array
c = SharedArray(eltype(rho), length(m));

simulate_dyn!(rho, c, vx, vy, vz, pos, m, A_FROM, A_TO)

save_realization(pos,vx,vy,vz,m,LOGGING_FOLDER)
