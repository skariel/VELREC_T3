@everywhere function _get_periodic_d(from, to, side_len=SIDE_LEN)
    d = to-from
    if d > side_len/2
        d -= side_len
    elseif d < -side_len/2
        d += side_len
    end
    d
end

@everywhere function _move_opos_i_inregards_to_pushed_target_single_worker!(opos_i, opos_f, _s_pos, frac_of_way, side_len=SIDE_LEN)
    for i in myrange(opos_i)
        d = _get_periodic_d(opos_f[i], _s_pos[i], side_len)
        opos_i[i] = mod1(opos_i[i] + frac_of_way * d, side_len)
    end
end

function _move_opos_i_inregards_to_pushed_target!(opos_i, opos_f, frac_of_way, side_len=SIDE_LEN)
    @sync begin
        for p in workers()
            @async remotecall_wait(p, _move_opos_i_inregards_to_pushed_target_single_worker!,
                opos_i, opos_f, _s_pos, frac_of_way, side_len)
        end
    end
end

function optimize_zeld_vs_pushed_pos!(rho, opos_i, pos, m, a_from, a_to, frac_mov=0.15, end_meandx=400.0, side_len=SIDE_LEN)
    info("optzel start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, " fracmov=",frac_mov)
    step = 0
    while true
        step += 1

        copy_into!(pos, opos_i)
        simulate_zeld!(rho, pos, m, a_from, a_to, side_len)

        (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
        info("optzel step=", step, " mdx=",mdx)
        mdx < end_meandx && break

        _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov, side_len)
    end
    info("optzel end")
end

function optimize_2lpt_vs_pushed_pos!(rho, opos_i, pos, m, a_from, a_to, frac_mov=0.15, end_meandx=400.0, side_len=SIDE_LEN)
    info("opt2lpt start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, " fracmov=",frac_mov)
    step = 0
    const fac2 = 1.0
    while true
        step += 1

        copy_into!(pos, opos_i)
        simulate_2lpt!(rho, pos, m, a_from, a_to, fac2, side_len)

        (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
        info("opt2lpt step=", step, " mdx=",mdx)
        mdx < end_meandx && break

        _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov, side_len)
    end
    info("opt2lpt end")
end

function optimize_dyn_vs_pushed_pos!(vx,vy,vz, rho, opos_i, pos, m, a_from, a_to, frac_mov=0.15, end_meandx=400.0)
    info("optdyn start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, " fracmov=",frac_mov)
    step = 0
    const fac2 = 1.0
    while true
        step += 1

        copy_into!(pos, opos_i)
        simulate_dyn!(rho, c, vx,vy,vz, pos, m, a_from, a_to)

        (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
        info("optdyn step=", step, " mdx=",mdx)
        mdx < end_meandx && break

        _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov, SIDE_LEN)
    end
    info("optdyn end")
end


function back_optimize_zeld_vs_pushed_pos!(rho, opos_i, pos, m, a_from, a_to, a_steps_num=20, frac_mov=0.15, end_meandx=400.0, side_len=SIDE_LEN)
    info("backoptzel start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, " fracmov=",frac_mov)
    lns = append!(collect(linspace(0.99,0.041,20)), (collect(linspace(0.04,0.02,120))))
    for a_i in lns
        if a_i > 0.3
            frac_mov = 0.27
        elseif a_i > 0.1
            frac_mov = 0.2
        elseif a_i > 0.05
            frac_mov = 0.1
        else a_i > 0.05
            frac_mov = 0.05
        end
        a_i==a_to && continue
        step = 0
        info("backoptzel a_i=",a_i)
        additional_frac_mov = 1.0
        while true
            step += 1

            copy_into!(pos, opos_i)
            simulate_zeld!(rho, pos, m, a_i, a_to, side_len)

            (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
            _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov*additional_frac_mov, side_len)
            additional_frac_mov *= 0.95
            info("backoptzel step=", step, " mdx=",mdx)
            mdx < end_meandx && break

        end
    end
    info("backoptzel end")
end

function back_optimize_2lpt_vs_pushed_pos!(rho, opos_i, pos, m, a_from, a_to, a_steps_num=20, frac_mov=0.15, end_meandx=400.0, side_len=SIDE_LEN)
    info("backopt2lpt start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, " fracmov=",frac_mov)
    lns = append!(collect(linspace(0.99,0.041,20)), (collect(linspace(0.04,0.02,120))))
    for a_i in lns
        if a_i > 0.3
            frac_mov = 0.27
        elseif a_i > 0.1
            frac_mov = 0.2
        elseif a_i > 0.05
            frac_mov = 0.1
        else a_i > 0.05
            frac_mov = 0.05
        end
        a_i==a_to && continue
        step = 0
        info("backopt2lpt a_i=",a_i)
        fac2 = 1.0
        additional_frac_mov = 1.0
        while true
            step += 1

            copy_into!(pos, opos_i)
            simulate_2lpt!(rho, pos, m, a_i, a_to, fac2, side_len)

            (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
            _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov*additional_frac_mov, side_len)
            additional_frac_mov *= 0.95
            info("backopt2lpt step=", step, " mdx=",mdx)
            mdx < end_meandx && break

        end
    end
    info("backopt2lpt end")
end

function back_optimize_dyn_vs_pushed_pos!(vx,vy,vz,rho, opos_i, pos, m, a_from, a_to, a_steps_num=20, frac_mov=0.15, end_meandx=400.0)
    info("backoptdyn start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, " fracmov=",frac_mov)
    lns = append!(collect(linspace(0.99,0.041,20)), (collect(linspace(0.04,0.02,120))))
    for a_i in lns
        if a_i > 0.3
            frac_mov = 0.27
        elseif a_i > 0.1
            frac_mov = 0.2
        elseif a_i > 0.05
            frac_mov = 0.1
        else a_i > 0.05
            frac_mov = 0.05
        end
        a_i==a_to && continue
        step = 0
        info("backoptdyn a_i=",a_i)
        fac2 = 1.0
        additional_frac_mov = 1.0
        while true
            step += 1

            copy_into!(pos, opos_i)
            simulate_dyn!(rho, c, vx,vy,vz, pos, m, a_i, a_to)

            (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
            _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov*additional_frac_mov, SIDE_LEN)
            additional_frac_mov *= 0.95
            info("backoptdyn step=", step, " mdx=",mdx)
            mdx < end_meandx && break

        end
    end
    info("backoptdyn end")
end
