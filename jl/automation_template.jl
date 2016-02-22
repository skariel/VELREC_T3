sim_name = "r"*REALIZATION_TYPE*
           "_p"*string(PROCS_NUM)*
           "_mr"*string(MASS_REAL)*
           "_opt"*OPT_TYPE*
           "_sim"*SIM_TYPE*
           "_ics"*ICS

jl_sim_name = sim_name*".jl"
run_sim_name = "RUN_"*sim_name
queue_sim_name = "QUEUE_"*sim_name

open(jl_sim_name, "w") do f
    write(f,"""

    rm("$(jl_sim_name)")
    rm("$(run_sim_name)")
    rm("$(queue_sim_name)")

    const REALIZATION_TYPE = "$(REALIZATION_TYPE)"
    const PROCS_NUM        = $(PROCS_NUM)
    const MASS_REAL        = $(MASS_REAL)
    const OPT_TYPE         = "$(OPT_TYPE)"
    const SIM_TYPE         = "$(SIM_TYPE)"
    const ICS              = "$(ICS)"
    const LOGGING_FOLDER   = "$(sim_name)"

    include("jl/init_by_params.jl")
    """)
end

open(run_sim_name,"w") do f
    write(f, """
    #!/bin/sh
    stdbuf -o0 -e0 julia/bin/julia $(jl_sim_name)
    """)
end

open(queue_sim_name,"w") do f
    write(f, """
    sbatch --partition=$(PARTITION) -N1 -n$(PROCS_NUM) --output=realization_$(sim_name)/log.out ./$(run_sim_name)
    """)
end

try
    rm( "realization_$(sim_name)", recursive=true)
catch e
end

try
    mkdir( "realization_$(sim_name)")
catch e
end

run(`chmod +x $(run_sim_name)`)
run(`chmod +x $(queue_sim_name)`)
run(`./$(queue_sim_name)`)
