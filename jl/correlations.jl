function _grid2reallinear(grid)
	const N = size(grid)[1]
    const N2 = div(N, 2)
	count = zeros(Int64, N)
	line = zeros(Complex128, N)
    for z in 0:N-1, y in 0:N-1, x in 0:N-1
        const kx = x<=N2 ? x : x - N
        const ky = y<=N2 ? y : y - N
        const kz = z<=N2 ? z : z - N
        k2 = kx*kx+ky*ky+kz*kz
		ix = round(Int, floor(sqrt(k2)))+1
		count[ix] += 1
		line[ix] += abs(grid[x+1,y+1,z+1])
	end
	real(line ./ count)
end

function get_correlation(grid, divide_by_n3=true)
    g = fft(grid)
    for i in eachindex(g)
        g[i] = abs2(g[i])
    end
    ifft!(g)
    if divide_by_n3
        return _grid2reallinear(g) / length(grid)
    end
    _grid2reallinear(g)
end
