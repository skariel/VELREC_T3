function _in_place_add!(a,v,dim)
    @inbounds for i in 1:size(a)[2]
        a[dim,i] += v
    end
    a
end

function in_place_multiply!(m, fac)
    @inbounds for i in eachindex(m)
        m[i] *= fac
    end
    m
end

function filter_mass_above(pos,vx,vy,vz,m,mf)
    N = 0
    @inbounds for mi in m
        mi < mf && continue
        N += 1
    end
    npos = SharedArray(eltype(pos), (3, N))
    nvx = SharedArray(eltype(vx), N)
    nvy = SharedArray(eltype(vy), N)
    nvz = SharedArray(eltype(vz), N)
    nm = SharedArray(eltype(m), N)
    ix = 0
    @inbounds for i in eachindex(m)
        m[i] < mf && continue
        ix += 1
        npos.s[:,ix] = pos[:,i]
        nvx.s[ix] = vx[i]
        nvy.s[ix] = vy[i]
        nvz.s[ix] = vz[i]
        nm.s[ix] = m[i]
    end
    npos,nvx,nvy,nvz,nm
end

function back_in_box!(pos, box_min=0.0, side_len=SIDE_LEN)
    @inbounds for i in eachindex(pos)
        pos[i] = mod1(pos[i], SIDE_LEN)
    end
end

function rho_to_1st_order_vel_pot!(rho)
    to_delta!(rho);
    to_g_fft!(rho);
    in_place_multiply!(rho, -1.0/4π);
    rho
end

function get_1st_order_vel!(c, a, dim, first_order_vel_pot)
    from_cic_dim2!(c,pos,first_order_vel_pot,dim);
    const fac = -D(a)*F(a)*Ha(a)
    @inbounds for i in eachindex(c)
        c[i] = fac * real(c[i])
    end
end

function firsrt_order_vel_pot_to_sencond_order!(vpot)
    to_tlpt_delta!(vpot);
    to_g_fft!(vpot);
    in_place_multiply!(vpot, -1.0/4π);
    vpot
end

function get_2nd_order_vel!(c, a, dim, second_order_vel_pot)
    from_cic_dim2!(c,pos,second_order_vel_pot,dim);
    const fac = D2(a)*F2(a)*Ha(a)
    @inbounds for i in eachindex(c)
        c[i] = fac * real(c[i])
    end
end
