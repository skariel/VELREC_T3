sim_name = "relname_"*RELNAME*"_"*
           "procs_"*string(PROCS_NUM)*"_"*
           "afrom_"*string(A_FROM)*"_"*
           "ato_"*string(A_TO)

jl_sim_name = sim_name*".jl"
run_sim_name = "RUN_"*sim_name
queue_sim_name = "QUEUE_"*sim_name

open(jl_sim_name, "w") do f
    write(f,"""

    rm("$(jl_sim_name)")
    rm("$(run_sim_name)")
    rm("$(queue_sim_name)")

    const A_FROM  = $(A_FROM)
    const A_TO    = $(A_TO)
    const RELNAME = "$(RELNAME)"
    const LOGGING_FOLDER   = "$(sim_name)"

    include("jl/init_by_params_sim.jl")
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
