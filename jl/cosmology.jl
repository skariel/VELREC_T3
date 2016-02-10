# using Gadget2 standard units:
#     velocity is canonical momentum, in km/s
#     distance is in kpc/h
#     mass is in 1e10 Ms/h
const G = 43007.1
const H0 = 0.1

# using Dark Sky cosmology:
const ΩΛ = 0.7048
const Ω0 = 1 - ΩΛ
const h = 0.6881
const σ8 = 0.8344

# using Dark-Sky halos box 1600 Mpc/h:
const SIDE_LEN = 1.6e6
const RHO_CRIT = 3*H0^2/8/π/G
const BOX_VOLUME = SIDE_LEN^3
const MASS_IN_BOX = RHO_CRIT*BOX_VOLUME*Ω0

# Hubble as function of scale factor a
Ha(a) = H0*sqrt(Ω0/a/a/a+ΩΛ)

# drift factor
FD(a1,a2) = quadgk(a->1./Ha(a)/a/a/a,a1,a2)[1]

# kick factor
FK(a1,a2) = quadgk(a->1./Ha(a)/a/a,a1,a2)[1]

# for convenience...
const FAC1 = 0.5*H0*H0*Ω0

# see here: http://mnras.oxfordjournals.org/content/322/2/419.full
OM(a) = Ω0 / (Ω0 + ΩΛ*a*a*a )
OL(a) = ΩΛ / (Ω0/a/a/a + ΩΛ )

function unnormalizedD(a)
    o0 = OM(a)
    ol = OL(a)
    up = 5o0
    down = o0^(4/7) - ol + (1+o0/2)*(1+ol/70)
    a*up/2down
end
D(a) = unnormalizedD(a) / unnormalizedD(1.0)
D2(a) = -3*D(a)^2/7

function F(a)
    o0 = OM(a)
    ol = OL(a)
    o0^(4/7) + (1+o0/2)*ol/70
end

F2(a) = 2*OM(a)^(6/11)

BOX_N = 512
SMTH = 3000.0
MEAN_SEP = 3010.314648627409
PARTICLE_NUM = 150149632
