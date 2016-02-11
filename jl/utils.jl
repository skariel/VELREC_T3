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
