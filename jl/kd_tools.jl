function get_slope_std(v_rec, v_real)
    const K = std(v_real) / std(v_rec)
    a,b = linreg(v_rec*K, v_real)
    b *= K
    b, std(v_real - v_rec*b)
end

function get_slope_std_smoothed(v_rec, v_real, pos, smth, n)
    ixs = rand(1:size(pos)[2], n);
    rngs = inrange(kd, pos[:,ixs], smth);
    _v_real_smth = zeros(length(rngs))
    _v_rec_smth = zeros(length(rngs))
    for i in eachindex(rngs)
        _v_real_smth[i] = mean(v_real[rngs[i]])
        _v_rec_smth[i] = mean(v_rec[rngs[i]])
    end
    get_slope_std(_v_rec_smth, _v_real_smth)
end

get_kd(pos) = KDTree(pos.s, reorder=false)

function measure_mass_density_th(kd, m, where_pos, smth=8000.0)
    const vol = 4/3*π*smth*smth*smth
    ixs = inrange(kd, where_pos, smth);
    sum(m[ixs])/vol
end

function measure_point_density_th(kd, where_pos, smth=8000.0)
    const vol = 4/3*π*smth*smth*smth
    ixs = inrange(kd, where_pos, smth);
    length(ixs)/vol
end

function get_random_point_in_box(box_min=0.0, side_len=SIDE_LEN, smth=8000.0)
    smth+rand(3)*(side_len-2*smth)
end
