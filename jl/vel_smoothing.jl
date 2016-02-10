# TODO: implement!


function get_slope_std(v_rec, v_real)
    const K = std(v_real) / std(v_rec)
    a,b = linreg(v_rec*K, v_real)
    b *= K
    b, std(v_real - v_rec*b)
end
