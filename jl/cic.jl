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

function to_cic!(x_arr,y_arr,z_arr,v_arr,grid,grid_min::Number,side_len::Number)
    fill!(grid, 0)
    @sync begin
        for p in workers()
            @async remotecall_wait(p, to_cic_single_worker!,
                x_arr, y_arr, z_arr, v_arr, grid, grid_min, side_len)
        end
    end
end

@everywhere function to_cic_single_worker!(x_arr,y_arr,z_arr,v_arr,grid,grid_min::Number,side_len::Number)
    fill!(grid, 0.0)
    const N = size(grid)[1]
    const g_dx = side_len / N
    @inbounds for i in myrange(x_arr)
        const x = mod1(x_arr[i] - grid_min, side_len)
        const y = mod1(y_arr[i] - grid_min, side_len)
        const z = mod1(z_arr[i] - grid_min, side_len)
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

function from_cic!(x_arr,y_arr,z_arr,grid,grid_min::Number,side_len::Number)
    v_arr = SharedArray(eltype(grid), length(x_arr))
    @sync begin
        for p in workers()
            @async remotecall_wait(p, from_cic_single_worker!,
                v_arr, x_arr, y_arr, z_arr, grid, grid_min, side_len)
        end
    end
    v_arr
end

@everywhere function from_cic_single_worker!(v_arr,x_arr,y_arr,z_arr,grid,grid_min::Number,side_len::Number)
    const N = size(grid)[1]
    const g_dx = side_len / eltype(x_arr)(N)
    @inbounds for i in myrange(v_arr)
        const x = mod1(x_arr[i] - grid_min, side_len)
        const y = mod1(y_arr[i] - grid_min, side_len)
        const z = mod1(z_arr[i] - grid_min, side_len)
        x_ix = 1 + round(Int, trunc(x/g_dx))
        y_ix = 1 + round(Int, trunc(y/g_dx))
        z_ix = 1 + round(Int, trunc(z/g_dx))
        v_arr[i] = 0
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

function in_place_add!(a,v)
    @inbounds for i in eachindex(a)
        a[i] += v
    end
    a
end

function from_cic_dx(x_arr,y_arr,z_arr,grid,grid_min::Number,side_len::Number)
    const dx = side_len / size(grid)[1]
    in_place_add!(x_arr, dx)   # x_arr + dx
    const pos1 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(x_arr, -2dx) # x_arr-dx
    const neg1 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(x_arr, 3dx) # x_arr+2dx
    const pos2 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(x_arr, -4dx) # x_arr-2dx
    const neg2 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    idx = 1.0/dx
    @inbounds for i in eachindex(pos1)
        pos1[i] = idx*(2/3*(pos1[i]-neg1[i])-1/12*(pos2[i]-neg2[i]))
    end
end

function from_cic_dy(x_arr,y_arr,z_arr,grid,grid_min::Number,side_len::Number)
    const dx = side_len / size(grid)[1]
    in_place_add!(y_arr, dx)   # x_arr + dx
    const pos1 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(y_arr, -2dx) # x_arr-dx
    const neg1 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(y_arr, 3dx) # x_arr+2dx
    const pos2 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(y_arr, -4dx) # x_arr-2dx
    const neg2 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    idx = 1.0/dx
    @inbounds for i in eachindex(pos1)
        pos1[i] = idx*(2/3*(pos1[i]-neg1[i])-1/12*(pos2[i]-neg2[i]))
    end
end

function from_cic_dz(x_arr,y_arr,z_arr,grid,grid_min::Number,side_len::Number)
    const dx = side_len / size(grid)[1]
    in_place_add!(z_arr, dx)   # x_arr + dx
    const pos1 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(z_arr, -2dx) # x_arr-dx
    const neg1 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(z_arr, 3dx) # x_arr+2dx
    const pos2 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    in_place_add!(z_arr, -4dx) # x_arr-2dx
    const neg2 = from_cic!(x_arr,y_arr,z_arr,grid,grid_min,side_len)
    idx = 1.0/dx
    @inbounds for i in eachindex(pos1)
        pos1[i] = idx*(2/3*(pos1[i]-neg1[i])-1/12*(pos2[i]-neg2[i]))
    end
end
