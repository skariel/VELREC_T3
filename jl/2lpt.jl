function to_tlpt_delta!(grid, side_len=SIDE_LEN)
    dx = similar(grid)
    dy = similar(grid)
    dz = similar(grid)
    const N = size(grid)[1]
    const GRID_DX = side_len/N
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        const nx = x+1 > N ? 1 : x+1
        const px = x-1 < 1 ? N : x-1
        dx[x,y,z] = (grid[nx,y,z] - grid[px,y,z])/2/GRID_DX
        const ny = y+1 > N ? 1 : y+1
        const py = y-1 < 1 ? N : y-1
        dy[x,y,z] = (grid[x,ny,z] - grid[x,py,z])/2/GRID_DX
        const nz = z+1 > N ? 1 : z+1
        const pz = z-1 < 1 ? N : z-1
        dz[x,y,z] = (grid[x,y,nz] - grid[x,y,pz])/2/GRID_DX
    end
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        const nx = x+1 > N ? 1 : x+1
        const px = x-1 < 1 ? N : x-1
        const ny = y+1 > N ? 1 : y+1
        const py = y-1 < 1 ? N : y-1
        const nz = z+1 > N ? 1 : z+1
        const pz = z-1 < 1 ? N : z-1

        const dxx = (dx[nx,y,z] - dx[px,y,z])/2/GRID_DX
        const dyy = (dy[x,ny,z] - dy[x,py,z])/2/GRID_DX
        const dzz = (dz[x,y,nz] - dz[x,y,pz])/2/GRID_DX

        const dyx = (dy[nx,y,z] - dy[px,y,z])/2/GRID_DX
        const dzx = (dz[nx,y,z] - dz[px,y,z])/2/GRID_DX
        const dzy = (dz[x,ny,z] - dz[x,py,z])/2/GRID_DX

        grid[x,y,z] = dyy*dxx + dzz*dxx + dzz*dyy  -  dyx*dyx - dzx*dzx - dzy*dzy
    end
end
