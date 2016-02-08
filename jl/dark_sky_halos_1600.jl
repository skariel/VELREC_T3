const RAW_DS_FN = "/home/skariel/Downloads/ds14_g_1600_4096_halos_1.0000"

immutable Halo1600
    x::Float32
    y::Float32
    z::Float32
    vx::Float32
    vy::Float32
    vz::Float32
    mvir::Float32
    m200b::Float32
    m200c::Float32
    m500c::Float32
    m2500c::Float32
    vmax::Float32
    rvmax::Float32
    r200b::Float32
    spin::Float32
    kin_to_pot::Float32
    id::Int64
    pid::Int64
end

function load_raw_halos1600(minimalmass=20.0, substructure=false)
    # translate mass to DS units...
    minimalmass *= 1e10

    const N = 3000000
    const HEADER_SIZE = 1216 + (20+4)*64
    tmp_halos = Array(Halo1600, N)
    chalos = Halo1600[]
    const sz = stat(RAW_DS_FN).size - HEADER_SIZE
    open(RAW_DS_FN, "r") do f
        # skip the header and sha1...
        seek(f, HEADER_SIZE)

        # eat large chunks of the file
        while position(f) < sz-sizeof(tmp_halos)
            # read halos
            read!(f, tmp_halos)

            # filter
            tot = 0
            for h in tmp_halos
                !substructure && h.pid >= 0 && continue
                h.m200b < minimalmass && continue
                tot += 1
                push!(chalos, h)
            end
            println("percentage read: ",(position(f)-HEADER_SIZE)/sz*100.0," fraction surviving mass filter: ",tot/N," halos#: ",length(chalos))
        end

        # eat whatever remains
        tmp_halos2 = Array(Halo1600, round(Int, (sz-position(f))/sizeof(Halo1600)))
        # read halos
        read!(f, tmp_halos2)

        if length(tmp_halos2) > 0
            # filter by substructure
            for h in tmp_halos2
                !substructure && h.pid >= 0 && continue
                h.m200b < minimalmass && continue
                push!(chalos, h)
            end
            println(position(f)/sz*100.0)
        end
    end
    # translate units to G2
    @inbounds for i in eachindex(chalos)
        chalos[i] = Halo1600(
            chalos[i].x*1000.0,
            chalos[i].y*1000.0,
            chalos[i].z*1000.0,
            chalos[i].vx,
            chalos[i].vy,
            chalos[i].vz,
            chalos[i].mvir/1e10,
            chalos[i].m200b/1e10,
            chalos[i].m200c/1e10,
            chalos[i].m500c/1e10,
            chalos[i].m2500c/1e10,
            chalos[i].vmax,
            chalos[i].rvmax,
            chalos[i].r200b*1000.0,
            chalos[i].spin,
            chalos[i].kin_to_pot,
            chalos[i].id,
            chalos[i].pid
        )
    end
    chalos
end

function get_pos_vel_mass(halos::Array{Halo1600,1})
    x = zeros(Float32, length(halos))
    y = zeros(Float32, length(halos))
    z = zeros(Float32, length(halos))
    vx = zeros(Float32, length(halos))
    vy = zeros(Float32, length(halos))
    vz = zeros(Float32, length(halos))
    m = zeros(Float32, length(halos))
    @inbounds for i in eachindex(halos)
        x[i] = halos[i].x
        y[i] = halos[i].y
        z[i] = halos[i].z
        vx[i] = halos[i].vx
        vy[i] = halos[i].vy
        vz[i] = halos[i].vz
        m[i] = halos[i].m200b
    end
    x,y,z,vx,vy,vz,m
end
