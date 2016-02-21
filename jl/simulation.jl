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
    TODO: FIX as above!!!
    info("sim2lpt start from a=",a_from," to a=",a_to)
    to_rho!(pos,m, rho);
    rho_to_1st_order_vel_pot!(rho);
    fac1 = 1.0
    for dim in 1:3
        info("sim2lpt 1st order dim ",dim)
        get_1st_order_s!(c, a_from, a_to, dim, rho)
        move_periodic!(pos, dim, c, fac1, side_len)
    end
    first_order_vel_pot_to_sencond_order!(rho)
    for dim in 1:3
        info("sim2lpt 2nd order dim ",dim)
        get_2nd_order_s!(c, a_from, a_to, dim, rho)
        move_periodic!(pos, dim, c, fac2, side_len)
    end
    info("sim2lpt end")
end

# doc"""
# returns dynamical allowed `da`
# this uses a few conservative approximations for speed:
# 1. 3D accel is approximated using a 1D accel
# 2. particles are skipped with the `skp` parameter
# """
# function get_min_dt(accel_x, a, smth=SMTH, skp=16)
#     minda = 1.e30 # infinity, ha!
#     hat = Ha(a)*a
#     const DTFRAC = 0.07
#     a4 = a*a*a*a
#     phys_smth = smth*a
#     @inbounds for i in 1:skp:length(accel_x)
#         const accel2 = accel_x[i]*accel_x[i]*3
#         const phys_accel2 = accel2/a4
#
#         const dyn_phys_dt2 = 2.0*DTFRAC*phys_smth/sqrt(phys_accel2)
#
#         const da = sqrt(dyn_phys_dt2)*hat
#
#         if da < minda
#             minda = da
#         end
#     end
#     if minda > 0.05
#         minda = 0.05
#     end
#     minda
# end
#
# function _kick_single_worker!(accel, v, a, da)
#     const fac1 = FK(a,a+da)
#     @inbounds for i in myrange(v)
#         v[i] += accel[i]*fac1
#     end
# end
#
# function kick!(accel, v, a, da)
#     @sync begin
#         for p in workers()
#             @async remotecall_wait(p, _kick_single_worker!,
#             accel, v, a, da)
#         end
#     end
# end
#
# function _drift_single_worker!(vel, pos, a, da, dim, side_len=SIDE_LEN)
#     const fac1 = FD(a,a+da)
#     @inbounds for i in myrange(v)
#         pos[dim,i] = mod1(pos[dim,i] + vel[i]*fac1, side_len)
#         pos[1,i] = mod1(pos[1,i] + real(vx[i])*fac, side_len)
#     end
# end
#
# function drift!(vel, pos, a, da, dim, side_len=SIDE_LEN)
#     @sync begin
#         for p in workers()
#             @async remotecall_wait(p, _drift_single_worker!,
#             vel, pos, a, da, dim, side_len)
#         end
#     end
# end
#
# function simulate_dyn!(rho, c, vx,vy,vz, pos, m, a_from, a_to, side_len=SIDE_LEN)
#     info("simdyn start from a=",a_from," to a=",a_to)
#     to_rho!(pos,m, rho);
#     rho_to_1st_order_vel_pot!(rho);
#     fac1 = 1.0
#     for dim in 1:3
#         info("simdyn dim ",dim)
#         get_1st_order_s!(c, a_from, a_to, dim, rho)
#         move_periodic!(pos, dim, c, fac1, side_len)
#     end
#     info("simdyn end")
# end
