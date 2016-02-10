FFTW.set_num_threads(8)

doc"""
Given a 3D cubed box represented by $\rho$, solves the poisson equation
$\nabla^2\phi=-4\pi\rho$ in a periodic fashion
"""
function to_g_fft!(ρ, side_len=SIDE_LEN, smth_len=SMTH)
    fft!(ρ)
    const N = size(ρ)[1]
    const N1 = N-1
    const N2 = div(N, 2)
    const fac = -N*N/π
    const fac2 = side_len/N
    asmth2 = 2*pi*smth_len/side_len;
    asmth2 *= asmth2;
    @inbounds for z in 0:N1, y in 0:N1, x in 0:N1
        const kx = x<=N2 ? x : x - N
        const ky = y<=N2 ? y : y - N
        const kz = z<=N2 ? z : z - N
        k2 = kx*kx+ky*ky+kz*kz
        if k2 > 0
            smth = -exp(-k2*asmth2); # TODO: remove this minus sign!
            ρ[x+1,y+1,z+1] *= fac/k2/fac2*smth
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
    @inbounds for i in CartesianRange(size(ρ))
        ρ[i] = ρ[i]/mn - 1.0
    end
    ρ
end

function in_place_multiply!(m, fac)
    @inbounds for i in eachindex(m)
        m[i] *= fac
    end
    m
end
