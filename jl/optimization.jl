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
    info("optzel start a_from=",a_from," a_to=",a_to," end_meandx=",end_meandx, "fracmov=",frac_mov)
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
