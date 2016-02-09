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

function save_realization(x,y,z,vx,vy,vz,m)
    # breaking data into chunks so github doesn't complain of files larger than 100Mb
    # github doesn't like it anyway... ignoring them in git now
    const N = length(x)
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

        save_arr(x[rng], "realization/x"*snum)
        save_arr(y[rng], "realization/y"*snum)
        save_arr(z[rng], "realization/z"*snum)
        save_arr(vx[rng], "realization/vx"*snum)
        save_arr(vy[rng], "realization/vy"*snum)
        save_arr(vz[rng], "realization/vz"*snum)
        save_arr(m[rng], "realization/m"*snum)

        i = f+1
    end
end

function load_chunked_array(fn)
    x = Float32[]
    for fix in 1:10
        snum = string(fix)*".gitignore"
        append!(x , load_arr_f32("realization/"*fn*snum))
    end
    sx = SharedArray(Float32, length(x))
    sx[:] = x
    sx
end

function load_realization()
    load_chunked_array("x"),
    load_chunked_array("y"),
    load_chunked_array("z"),
    load_chunked_array("vx"),
    load_chunked_array("vy"),
    load_chunked_array("vz"),
    load_chunked_array("m")
end
