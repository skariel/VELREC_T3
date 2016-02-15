
function to_tlpt_delta!(grid, side_len=SIDE_LEN)
    res = similar(grid)
    fill!(res, 0.0)
    const N = size(grid)[1]
    const GRID_DX = side_len/N
    const fac = 16*GRID_DX^4

    @inbounds @fastmath for z in 1:N
        for y in 1:N
            @simd for x in 1:N

                const xp1 = mod1(x+1,N)
                const xm1 = mod1(x-1,N)
                const xp2 = mod1(x+2,N)
                const xm2 = mod1(x-2,N)

                const yp1 = mod1(y+1,N)
                const ym1 = mod1(y-1,N)
                const yp2 = mod1(y+2,N)
                const ym2 = mod1(y-2,N)

                const zp1 = mod1(z+1,N)
                const zm1 = mod1(z-1,N)
                const zp2 = mod1(z+2,N)
                const zm2 = mod1(z-2,N)

                const g02 = 2*grid[x,y,z]
                const dxx = grid[xp2,y,z] - g02 + grid[xm2,y,z]
                const dyy = grid[x,yp2,z] - g02 + grid[x,ym2,z]
                const dzz = grid[x,y,zp2] - g02 + grid[x,y,zm2]

                const dyx = grid[xp1,yp1,z]-grid[xp1,ym1,z] - (grid[xm1,yp1,z]-grid[xm1,ym1,z])
                const dzx = grid[xp1,y,zp1]-grid[xp1,y,zm1] - (grid[xm1,y,zp1]-grid[xm1,y,zm1])
                const dzy = grid[x,yp1,zp1]-grid[x,yp1,zm1] - (grid[x,ym1,zp1]-grid[x,ym1,zm1])

                res[x,y,z] = (dyy*dxx + dzz*dxx + dzz*dyy  -  dyx*dyx - dzx*dzx - dzy*dzy)/fac
            end
        end
    end

    @inbounds @fastmath for z in 1:N
        for y in 1:N
            @simd for x in 1:N
                grid[x,y,z] = res[x,y,z]
            end
        end
    end

end
