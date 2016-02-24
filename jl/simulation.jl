function simulate_zeld!(rho, pos, m, a_from, a_to, side_len=SIDE_LEN)
    info("simzel start from a=",a_from," to a=",a_to)
    to_rho!(pos,m, rho);
    rho_to_1st_order_vel_pot!(rho);
    fac1 = 1.0

    info("simzel dim 1")
    get_1st_order_s!(c, a_from, a_to, 1, pos, rho)
    dx = real(c)
    info("simzel dim 2")
    get_1st_order_s!(c, a_from, a_to, 2, pos, rho)
    dy = real(c)
    info("simzel dim 3")
    get_1st_order_s!(c, a_from, a_to, 3, pos, rho)
    dz = real(c)

    move_periodic!(pos, 1, dx, fac1, side_len)
    move_periodic!(pos, 2, dy, fac1, side_len)
    move_periodic!(pos, 3, dz, fac1, side_len)

    info("simzel end")
end

function simulate_2lpt!(rho, pos, m, a_from, a_to, fac2=1.0, side_len=SIDE_LEN)
    # TEST!
    info("sim2lpt start from a=",a_from," to a=",a_to)
    to_rho!(pos,m, rho);
    rho_to_1st_order_vel_pot!(rho);

    info("sim2lpt 1st order dim 1")
    get_1st_order_s!(c, a_from, a_to, 1, pos, rho)
    dx = real(c)
    info("sim2lpt 1st order dim 2")
    get_1st_order_s!(c, a_from, a_to, 2, pos, rho)
    dy = real(c)
    info("sim2lpt 1st order dim 3")
    get_1st_order_s!(c, a_from, a_to, 3, pos, rho)
    dz = real(c)

    first_order_vel_pot_to_sencond_order!(rho)

    info("sim2lpt 2nd order dim 1")
    get_2nd_order_s!(c, a_from, a_to, 1, pos, rho)
    dx += real(c)*fac2
    info("sim2lpt 2nd order dim 2")
    get_2nd_order_s!(c, a_from, a_to, 2, pos, rho)
    dy += real(c)*fac2
    info("sim2lpt 2nd order dim 3")
    get_2nd_order_s!(c, a_from, a_to, 3, pos, rho)
    dz += real(c)*fac2

    move_periodic!(pos, 1, dx, 1.0, side_len)
    move_periodic!(pos, 2, dy, 1.0, side_len)
    move_periodic!(pos, 3, dz, 1.0, side_len)
    info("sim2lpt end")
end

doc"""
returns dynamical allowed `da`
this uses a few conservative approximations for speed:
1. 3D accel is approximated using a 1D accel
2. particles are skipped with the `skp` parameter
"""
function get_min_dt(accel_x, a, smth=SMTH, skp=16)
    minda = 1.e30 # infinity, ha!
    hat = Ha(a)*a
    const DTFRAC = 0.07
    a4 = a*a*a*a
    phys_smth = smth*a
    das = Float64[]
    for i in 1:skp:length(accel_x)
        const accel2 = abs2(accel_x[i])*3
        const phys_accel2 = accel2/a4

        const dyn_phys_dt2 = 2.0*DTFRAC*phys_smth/sqrt(phys_accel2)

        const da = sqrt(dyn_phys_dt2)*hat

        # if da < minda
        #     minda = da
        # end
        push!(das, da)
    end
    minda = sort(das)[div(length(das), 20)]
    if minda > 0.2
        minda = 0.2
    end
    minda
end

@everywhere function _kick_single_worker!(accel, v, fk)
    for i in myrange(v)
        v[i] += real(accel[i])*fk
    end
end

function kick!(accel, v, a, da)
    const fk = FK(a,a+da)
    @sync begin
        for p in workers()
            @async remotecall_wait(p, _kick_single_worker!,
            accel, v, fk)
        end
    end
end

@everywhere function _drift_single_worker!(vel, pos, fd, dim, side_len=SIDE_LEN)
    for i in myrange(vel)
        pos[dim,i] = mod1(pos[dim,i] + vel[i]*fd, side_len)
    end
end

function drift!(vel, pos, a, da, dim, side_len=SIDE_LEN)
    const fd = FD(a,a+da)
    @sync begin
        for p in workers()
            @async remotecall_wait(p, _drift_single_worker!,
            vel, pos, fd, dim, side_len)
        end
    end
end

function simulate_dyn!(rho, c, vx,vy,vz, pos, m, a_from, a_to)
    info("simdyn start from a=",a_from," to a=",a_to)

    # Initial velocity
    info("simdyn initial vel")
    to_rho!(pos,m, rho);
    rho_to_1st_order_vel_pot!(rho);
    get_1st_order_comoving_vel!(c, a_from, 1, pos, rho)
    vx.s[:] = real(c)*a_from*a_from
    get_1st_order_comoving_vel!(c, a_from, 2, pos, rho)
    vy.s[:] = real(c)*a_from*a_from
    get_1st_order_comoving_vel!(c, a_from, 3, pos, rho)
    vz.s[:] = real(c)*a_from*a_from

    # Initial accel...
    info("simdyn initial accel")
    to_rho!(pos,m, rho)
    to_g_fft!(rho)
    in_place_multiply!(rho ,G)
    from_cic_dim2!(c,pos,rho,1)
    step = 0
    a = a_from

    finished = false
    while !finished

        # calc da
        da = get_min_dt(c,a,SMTH,16)
        if a+da > a_to
            finished = true
            da = a_to - a
        end

        info("simdyn step=",step," a=",a," da=",da)

        # kick da/2
        info("simdyn kick da/2.. x")
        from_cic_dim2!(c,pos,rho,1)
        kick!(c, vx, a, da/2)
        info("simdyn kick da/2.. y")
        from_cic_dim2!(c,pos,rho,2)
        kick!(c, vy, a, da/2)
        info("simdyn kick da/2.. z")
        from_cic_dim2!(c,pos,rho,3)
        kick!(c, vz, a, da/2)

        # drift da
        info("simdyn drift da")
        drift!(vx, pos, a, da, 1)
        drift!(vy, pos, a, da, 2)
        drift!(vz, pos, a, da, 3)

        # accel
        info("simdyn accel")
        to_rho!(pos,m, rho)
        to_g_fft!(rho)
        in_place_multiply!(rho ,G)

        # kick da/2
        info("simdyn kick da/2.. x")
        from_cic_dim2!(c,pos,rho,1)
        kick!(c, vx, a, da/2)
        info("simdyn kick da/2.. y")
        from_cic_dim2!(c,pos,rho,2)
        kick!(c, vy, a, da/2)
        info("simdyn kick da/2.. z")
        from_cic_dim2!(c,pos,rho,3)
        kick!(c, vz, a, da/2)

        # just advance time...
        a += da
        step += 1
    end

    info("simdyn end")
end
