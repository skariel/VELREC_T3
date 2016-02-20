# TODO: implement!

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
    @inbounds for i in myrange(opos_i)
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

        pos.s[:,:] = opos_i.s[:,:]
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

        pos.s[:,:] = opos_i.s[:,:]
        simulate_2lpt!(rho, pos, m, a_from, a_to, fac2, side_len)

        (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
        info("opt2lpt step=", step, " mdx=",mdx)
        mdx < end_meandx && break

        _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov, side_len)
    end
    info("opt2lpt end")
end

function back_optimize_zeld_vs_pushed_pos!(rho, opos_i, pos, m, a_from, a_to, a_steps_num=20, frac_mov=0.15, end_meandx=400.0, side_len=SIDE_LEN)
    info("backoptzel start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, " fracmov=",frac_mov)
    for a_i in linspace(a_to,a_from,a_steps_num)
        a_i==a_to && continue
        step = 0
        info("backoptzel a_i=",a_i)
        while true
            step += 1

            pos.s[:,:] = opos_i.s[:,:]
            simulate_zeld!(rho, pos, m, a_i, a_to, side_len)

            (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
            _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov, side_len)
            info("backoptzel step=", step, " mdx=",mdx)
            mdx < end_meandx && break

        end
    end
    info("backoptzel end")
end

function back_optimize_2lpt_vs_pushed_pos!(rho, opos_i, pos, m, a_from, a_to, a_steps_num=20, frac_mov=0.15, end_meandx=400.0, side_len=SIDE_LEN)
    info("backopt2lpt start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, " fracmov=",frac_mov)
    for a_i in linspace(a_to,a_from,a_steps_num)
        a_i==a_to && continue
        step = 0
        info("backopt2lpt a_i=",a_i)
        fac2 = 1.0
        while true
            step += 1

            pos.s[:,:] = opos_i.s[:,:]
            simulate_2lpt!(rho, pos, m, a_i, a_to, fac2, side_len)

            (mdx,sdx) = mean_std_dx_vs_pushed_pos(pos)
            _move_opos_i_inregards_to_pushed_target!(opos_i, pos, frac_mov, side_len)
            info("backopt2lpt step=", step, " mdx=",mdx)
            mdx < end_meandx && break

        end
    end
    info("backopt2lpt end")
end
