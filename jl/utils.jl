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

function move_periodic_all_dims!(pos, vx, vy, vz, fac, side_len=SIDE_LEN)
    @inbounds for i in 1:length(vx)
        pos[1,i] = mod1(pos[1,i] + real(vx[i])*fac, side_len)
        pos[2,i] = mod1(pos[2,i] + real(vy[i])*fac, side_len)
        pos[3,i] = mod1(pos[3,i] + real(vz[i])*fac, side_len)
    end
end

function mean_std_dx_vs_pushed_pos(pos, side_len=SIDE_LEN)
    dx = zeros(Float32, div(size(pos)[2], 50))
    @inbounds for i in 1:length(dx)
        i50 = i*50
        dx[i] = abs(pos[1,i50]-_s_pos[1,i50])
        if dx[i] > side_len/2
            dx[i] = abs(dx[i]-side_len)
        end
    end
    mean(dx), std(dx)
end

function move_periodic!(pos, dim, vel, fac, side_len=SIDE_LEN)
    @inbounds for i in 1:length(vel)
        pos[dim,i] = mod1(pos[dim,i] + real(vel[i])*fac, side_len)
    end
end

function rho_to_1st_order_vel_pot!(rho)
    to_delta!(rho);
    to_g_fft!(rho);
    in_place_multiply!(rho, -1.0/4π);
    rho
end

function get_1st_order_comoving_vel!(c, a, dim, first_order_vel_pot)
    from_cic_dim2!(c,pos,first_order_vel_pot,dim);
    const fac = -D(a)*F(a)*Ha(a)
    @inbounds for i in eachindex(c)
        c[i] = fac * real(c[i])
    end
end

function get_1st_order_s!(c, a_from, a_to, dim, first_order_vel_pot)
    from_cic_dim2!(c,pos,first_order_vel_pot,dim);
    const fac_from = -D(a_from)
    const fac_to = -D(a_to)
    const dfac = fac_to-fac_from
    @inbounds for i in eachindex(c)
        c[i] = dfac * real(c[i])
    end
end

function move_1st_order!(pos, c, a_from, a_to, first_order_vel_pot, fac=1.0, side_len=SIDE_LEN)
    for dim in 1:3
        get_1st_order_s!(c, a_from, a_to, 1, first_order_vel_pot)
        move_periodic!(pos, dim, c, fac, side_len)
    end
end

function first_order_vel_pot_to_sencond_order!(vpot)
    to_tlpt_delta!(vpot);
    to_g_fft!(vpot);
    in_place_multiply!(vpot, -1.0/4π);
    vpot
end

function get_2nd_order_comoving_vel!(c, a, dim, second_order_vel_pot)
    from_cic_dim2!(c,pos,second_order_vel_pot,dim);
    const fac = D2(a)*F2(a)*Ha(a)
    @inbounds for i in eachindex(c)
        c[i] = fac * real(c[i])
    end
end

function get_2nd_order_s!(c, a_from, a_to, dim, second_order_vel_pot)
    from_cic_dim2!(c,pos,second_order_vel_pot,dim);
    const fac_from = D2(a_from)
    const fac_to = D2(a_to)
    const dfac = fac_to-fac_from
    @inbounds for i in eachindex(c)
        c[i] = dfac * real(c[i])
    end
end

function move_2nd_order!(pos, c, a_from, a_to, second_order_vel_pot, fac=1.0, side_len=SIDE_LEN)
    for dim in 1:3
        get_2nd_order_s!(c, a_from, a_to, 1, second_order_vel_pot)
        move_periodic!(pos, dim, c, fac, side_len)
    end
end
