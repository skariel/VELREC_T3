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

function load_raw_halos1600(minimalmass=10.0, substructure=false)
    # translate mass to DS units...
    minimalmass *= 1e10

    const N = 3000000
    const HEADER_SIZE = 1216 + (20+4)*64
    tmp_halos = Array(Halo1600, N)
    _x = Float32[]
    _y = Float32[]
    _z = Float32[]
    _vx = Float32[]
    _vy = Float32[]
    _vz = Float32[]
    _m = Float32[]
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
                push!(_x, h.x)
                push!(_y, h.y)
                push!(_z, h.z)
                push!(_vx, h.vx)
                push!(_vy, h.vy)
                push!(_vz, h.vz)
                push!(_m, h.m200b)
            end
            info("percentage read: ",(position(f)-HEADER_SIZE)/sz*100.0," fraction surviving mass filter: ",tot/N," halos#: ",length(_x))
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
                push!(_x, h.x)
                push!(_y, h.y)
                push!(_z, h.z)
                push!(_vx, h.vx)
                push!(_vy, h.vy)
                push!(_vz, h.vz)
                push!(_m, h.m200b)
            end
            info(position(f)/sz*100.0)
        end
    end
    # translate units to G2
    _x *= 1000.0
    _y *= 1000.0
    _z *= 1000.0
    _m /= 1f10
    _pos = zeros(eltype(_x), (3, length(_x)))
    _pos[1,:] = _x
    _pos[2,:] = _y
    _pos[3,:] = _z
    _pos, _vx, _vy, _vz, _m
end
