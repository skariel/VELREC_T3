function myrange(q)
    1:length(q)
end

@everywhere function myrange(q::SharedArray)
    idx = indexpids(q)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0,length(q),nchunks+1)]
    splits[idx]+1:splits[idx+1]
end

function to_cic!(data,v_arr,grid,grid_min=0.0,side_len=SIDE_LEN)
    fill!(grid, 0)
    @sync begin
        for p in workers()
            @async remotecall_wait(p, _to_cic_single_worker!,
                data, v_arr, grid, grid_min, side_len)
        end
    end
end

function to_rho!(data,v_arr,grid,grid_min=0.0,side_len=SIDE_LEN)
    const N = size(grid)[1]
    const DX = side_len/N
    to_cic!(data,v_arr,grid,grid_min,side_len)
    in_place_multiply!(rho, 1.0/DX/DX/DX)
end

@everywhere function _to_cic_single_worker!(data,v_arr,grid,grid_min::Number,side_len::Number)
    const N = size(grid)[1]
    const g_dx = side_len / eltype(data)(N)
    @inbounds for i in myrange(v_arr)
        const x = mod1(data[1,i] - grid_min, side_len)
        const y = mod1(data[2,i] - grid_min, side_len)
        const z = mod1(data[3,i] - grid_min, side_len)
        const x_ix = 1 + round(Int, trunc(x/g_dx))
        const y_ix = 1 + round(Int, trunc(y/g_dx))
        const z_ix = 1 + round(Int, trunc(z/g_dx))
        val = v_arr[i]
        @inbounds for d_ix_x in 0:1, d_ix_y in 0:1, d_ix_z in 0:1
            const lever_x = abs(x_ix-d_ix_x - x/g_dx)
            const lever_y = abs(y_ix-d_ix_y - y/g_dx)
            const lever_z = abs(z_ix-d_ix_z - z/g_dx)
            const _x_ix = x_ix + d_ix_x > N ? 1 : x_ix + d_ix_x
            const _y_ix = y_ix + d_ix_y > N ? 1 : y_ix + d_ix_y
            const _z_ix = z_ix + d_ix_z > N ? 1 : z_ix + d_ix_z
            grid[_x_ix, _y_ix, _z_ix] += val*lever_x*lever_y*lever_z
        end
    end
end

function from_cic!(v_arr,data,grid,grid_min=0.0,side_len=SIDE_LEN)
    fill!(v_arr, 0)
    @sync begin
        for p in workers()
            @async remotecall_wait(p, _from_cic_single_worker!,
                v_arr, data, grid, grid_min, side_len)
        end
    end
    v_arr
end

@everywhere function _from_cic_single_worker!(v_arr,data,grid,grid_min::Number,side_len::Number)
    const N = size(grid)[1]
    const g_dx = side_len / eltype(data)(N)
    @inbounds for i in myrange(v_arr)
        const x = mod1(data[1,i] - grid_min, side_len)
        const y = mod1(data[2,i] - grid_min, side_len)
        const z = mod1(data[3,i] - grid_min, side_len)
        x_ix = 1 + round(Int, trunc(x/g_dx))
        y_ix = 1 + round(Int, trunc(y/g_dx))
        z_ix = 1 + round(Int, trunc(z/g_dx))
        @inbounds for d_ix_x in 0:1, d_ix_y in 0:1, d_ix_z in 0:1
            const lever_x = abs(x_ix-d_ix_x - x/g_dx)
            const lever_y = abs(y_ix-d_ix_y - y/g_dx)
            const lever_z = abs(z_ix-d_ix_z - z/g_dx)
            const _x_ix = x_ix + d_ix_x > N ? 1 : x_ix + d_ix_x
            const _y_ix = y_ix + d_ix_y > N ? 1 : y_ix + d_ix_y
            const _z_ix = z_ix + d_ix_z > N ? 1 : z_ix + d_ix_z
            v_arr[i] += grid[_x_ix, _y_ix, _z_ix] * lever_x*lever_y*lever_z
        end
    end
end

function from_cic_dim!(v_arr,data,grid,dim::Integer, grid_min=0.0, side_len=SIDE_LEN)
    const dx = side_len / size(grid)[1]
    const idx = 1.0/dx
    cum_tmp = zeros(eltype(grid), length(v_arr))

    _in_place_add!(data, dx, dim)   # +dx
    from_cic!(v_arr,data,grid,grid_min,side_len)
    @inbounds for i in eachindex(v_arr)
        cum_tmp[i] += v_arr[i]*idx*2/3
    end
    _in_place_add!(data, -2dx, dim) # -dx
    tmp = from_cic!(v_arr,data,grid,grid_min,side_len)
    @inbounds for i in eachindex(v_arr)
        cum_tmp[i] -= v_arr[i]*idx*2/3
    end
    _in_place_add!(data, 3dx, dim)  # +2dx
    tmp = from_cic!(v_arr,data,grid,grid_min,side_len)
    @inbounds for i in eachindex(v_arr)
        cum_tmp[i] += v_arr[i]*idx* ( -1/12)
    end
    _in_place_add!(data, -4dx, dim) # -2dx
    tmp = from_cic!(v_arr,data,grid,grid_min,side_len)
    @inbounds for i in eachindex(v_arr)
        cum_tmp[i] -= v_arr[i]*idx* (-1/12)
    end

    # finalizing...
    _in_place_add!(data, 2dx, dim) # back to original...
    @inbounds for i in eachindex(v_arr)
        v_arr[i] = cum_tmp[i]
    end
    v_arr
end
