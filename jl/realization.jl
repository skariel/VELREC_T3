function save_arr(arr, fn)
    open(fn, "w") do f
        write(f,arr)
    end
end

function load_arr_f32(fn)
    const SZ = stat(fn).size
    const N = div(SZ, sizeof(Float32))
    arr = zeros(Float32, N)
    open(fn, "r") do f
        read!(f,arr)
    end
    arr
end

function save_realization(pos,vx,vy,vz,m)
    # breaking data into chunks so github doesn't complain of files larger than 100Mb
    # github doesn't like it anyway... ignoring them in git now
    const N = length(m)
    i = 1
    f = 1
    dix = div(N, 10)
    fix = 1
    while f<N
        snum = string(fix)*".gitignore"
        fix += 1

        f = i+dix
        if (f>N)
            f = N # this will be the last iteration
        end
        rng = i:f

        save_arr(pos[1,rng], "realization/x"*snum)
        save_arr(pos[2,rng], "realization/y"*snum)
        save_arr(pos[3,rng], "realization/z"*snum)
        save_arr(vx[rng], "realization/vx"*snum)
        save_arr(vy[rng], "realization/vy"*snum)
        save_arr(vz[rng], "realization/vz"*snum)
        save_arr(m[rng], "realization/m"*snum)

        i = f+1
    end
end

function load_chunked_array(fn, shared=true)
    x = Float32[]
    for fix in 1:10
        snum = string(fix)*".gitignore"
        append!(x , load_arr_f32("realization/"*fn*snum))
    end
    if shared
        sx = SharedArray(Float32, length(x))
        sx[:] = x
        return sx
    else
        sx = Array(Float32, length(x))
        sx[:] = x
        return sx
    end
end

function load_realization()
    x = load_chunked_array("x", false)
    y = load_chunked_array("y", false)
    z = load_chunked_array("z", false)
    pos = SharedArray(Float32, (3,length(x)))
    pos[1,:] = x
    pos[2,:] = y
    pos[3,:] = z

    pos,
    load_chunked_array("vx"),
    load_chunked_array("vy"),
    load_chunked_array("vz"),
    load_chunked_array("m")
end
