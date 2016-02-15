FFTW.set_num_threads(8)

doc"""
Given a 3D cube box represented by $\rho$, solves the poisson equation
$\nabla^2\phi=-4\pi\rho$ in a periodic fashion. Note \rho here is actual
density!
"""
function to_g_fft!(ρ, side_len=SIDE_LEN, smth_len=SMTH)
    fft!(ρ)
    const N = size(ρ)[1]
    const N2 = div(N, 2)
    const fac = side_len*side_len/π
    const ksmth2 = (2π*smth_len/side_len)^2
    @fastmath @inbounds for z in 0:N-1, y in 0:N-1, x in 0:N-1
        const kx = x<=N2 ? x : x - N
        const ky = y<=N2 ? y : y - N
        const kz = z<=N2 ? z : z - N
        k2 = kx*kx+ky*ky+kz*kz
        if k2 > 0
            ρ[x+1,y+1,z+1] *= fac/k2*exp(-k2*ksmth2)
        end
    end
    ifft!(ρ)
    const mn = ρ[end,end,end]
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        ρ[x,y,z] -= mn
    end
    ρ
end

function to_delta!(ρ)
    mn = mean(ρ)
    @fastmath @inbounds for i in CartesianRange(size(ρ))
        ρ[i] = ρ[i]/mn - 1.0
    end
    ρ
end
