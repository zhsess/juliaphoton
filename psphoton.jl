# Program to model seismic phases as energy packets
# Originally written 2003-2004 by Peter Shearer (pshearer@ucsd.edu)
# Converted to f90 in 2012 by Nick Mancinelli (n.j.mancinelli@gmail.com)
# Converted to Julia in 2024 by Hao Zhang (zhsess@gmail.com)
# Last modified: 2024-04-24

using SpecialFunctions
include("gavg.jl")
include("flatten.jl")
include("layertrace.jl")

# read in parameters
vmodel::String = "./iasp91.smo"
zsource::Float64 = 0.0 # assumed source depth (km)
freq::Float64 = 1.0 # assumed freq. (Hz), affects attenuation and scattering
om::Float64 = 2.0 * pi * freq

nqlay::Int64 = 3 # number of Q layers (max=6)
zminq = Float64[0, 220, 5153.9] # min depth of Q layer
zmaxq = Float64[220, 2889, 9999.0] # max depth of Q layer
qalpha = Float64[227, 1383, 360] # Q alpha

nscatvol::Int64 = 4 # number of scattering volumes (to follow) (max=6)
zminscat = Float64[1, 200, 600, 1700] # min scat depth (km)
zmaxscat = Float64[200, 600, 1700, 2880] # max scat depth (km)
xmaxscat = Float64[999999, 999999, 999999, 999999] # max scat range from source (km)
pvelref = Float64[8.08, 9.13, 11.734, 13.131] # reference P & S velocity for layer
svelref = Float64[4.47, 4.93, 6.504, 7.055]
el = om ./ svelref
gam0 = pvelref ./ svelref
nu = Float64[0.8, 0.8, 0.8, 0.8] # relative size of density perturbation (0.8 often assumed)
eps = Float64[0.02, 0.02, 0.005, 0.005] # rms perturbation
alen = Float64[4.0, 4.0, 8.0, 8.0] # scale length (km)

nface::Int64 = 3 # number of interface depths (to follow) (max=6)
depface = Float64[0, 2889, 5153.9] # iasp91 CMB & ICB

np::Int64 = 10000 # number of ray parameters for tables (max=50000)

# define variables
nface0::Int64 = 6
nlay0::Int64 = 655
nray0::Int64 = 50000
erad::Float64 = 6371.0

z_s = zeros(Float64, nlay0)
alpha_s = zeros(Float64, nlay0)
beta_s = zeros(Float64, nlay0)
rho = zeros(Float64, nlay0)

z = zeros(Float64, nlay0)
alpha = zeros(Float64, nlay0)
beta = zeros(Float64, nlay0)

q = zeros(Float64, nlay0, 2)
slow = zeros(Float64, nlay0, 2)

iflag = zeros(Int64, nlay0)
iflagvol = zeros(Int64, nlay0)

dep::Float64 = 0.0
izsource::Int64 = 0
npts::Int64 = 0

scatprob = zeros(Float64, nface0, 2, 3)

vp = zeros(Float64, nface, 2)
vs = zeros(Float64, nface, 2)
den = zeros(Float64, nface, 2)

amp0 = zeros(Float64, nray0, 2)
p = zeros(Float64, nray0, 2)
iddir = zeros(Int64, nlay0, nray0, 2)
dx = zeros(Float64, nlay0, nray0, 2)
dt = zeros(Float64, nlay0, nray0, 2)
dtstar = zeros(Float64, nlay0, nray0, 2)

#---------------------------------------------------------
# read velocity model and transform to flat earth
# value at center of earth is removed to avoid singularity

open(vmodel, "r") do f
    for i in 1:nlay0
        # Read the line and split it into variables
        line_ = split(readline(f))
        z_s[i] = parse(Float64, line_[1])
        dep = parse(Float64, line_[2])
        alpha_s[i] = parse(Float64, line_[3])
        beta_s[i] = parse(Float64, line_[4])
        rho[i] = parse(Float64, line_[5])

        if z_s[i] == zsource global izsource = i end
        if z_s[i] == erad
            global npts = i
            break
        end
        if beta_s[i] == 0.0 beta_s[i] = 0.001 end

        z[i], alpha[i] = FLATTEN(z_s[i], alpha_s[i])
        if beta_s[i] != 0.0
            z[i], beta[i] = FLATTEN(z_s[i], beta_s[i])
        else
            beta[i] = 0.0
        end
        iflag[i] = 0
        iflagvol[i] = 0

        q[i, 1] = 999999.
        for k in 1:nqlay
            if (z_s[i] >= zminq[k]) & (z_s[i] <= zmaxq[k])
                q[i, 1] = qalpha[k]
            end
        end
        q[i, 2] = (4/9) * q[i, 1]
        global npts = i
    end
end

vpmin = alpha[1]
vsmin = beta[1]
vpsource = alpha[izsource]
vssource = beta[izsource]

z_s[npts] = z_s[npts - 1]
alpha_s[npts] = alpha_s[npts - 1]
beta_s[npts] = beta_s[npts - 1]
q[npts, 1] = q[npts - 1, 1]
q[npts, 2] = q[npts - 1, 2]
z[npts], alpha[npts] = FLATTEN(z_s[npts], alpha_s[npts])
z[npts], beta[npts] = FLATTEN(z_s[npts], beta_s[npts])

println("finished reading model")
println("Depth points in model = $npts \n")

for i in 1:npts
    slow[i, 1] = 1.0 / alpha[i]
    if beta[i] != 0.0
        slow[i, 2] = 1.0 / beta[i]
    else
        slow[i, 2] = 1.0 / alpha[i]
    end
end

# set up scattering volume flags and other stuff
# Note that if scattering volumes overlap, iflagvol is set to last of scattering volumes
for iscat in 1:nscatvol
    gpp0, gps0, gsp0, gss0, gppm, gpsm, gspm, gssm = GAVGSATO2(el[iscat], nu[iscat], gam0[iscat], eps[iscat], alen[iscat])

    # println("Scattering Strengths for layer $iscat:")
    # println("gpp0 = $gpp0, gps0 = $gps0, gsp0 = $gsp0, gss0 = $gss0")
    # println("afp = ", 1. / (gpp0 + gps0), ", ", 1. / (gsp0 + gss0))
    # println("gppm = $gppm, gpsm = $gpsm, gspm = $gspm, gssm = $gssm")
    # println("moment afp =", 1. / (gppm + gpsm), ", ", 1. / (gspm + gssm), "\n")

    scatprob[iscat, 1, 1] = gpp0 + gps0
    scatprob[iscat, 1, 2] = gpp0
    scatprob[iscat, 1, 3] = gps0
    scatprob[iscat, 2, 1] = gsp0 + gss0
    scatprob[iscat, 2, 2] = gsp0
    scatprob[iscat, 2, 3] = gss0

    for i in 2:npts
        if (z_s[i-1] >= zminscat[iscat]) & (z_s[i] <= zmaxscat[iscat])
            iflagvol[i] = iscat
        end
    end
end

# -------------------- Now define interface values --------------------
depface[1] = 0.0
vp[1, 1] = 0.1
vs[1, 1] = 0.01
den[1, 1] = 0.01
vp[1, 2] = alpha_s[1]
vs[1, 2] = beta_s[1]
den[1, 2] = rho[1]
iflag[1] = 1

for iface in 2:nface
    for i = 2:npts
        if (z_s[i] == z_s[i-1]) & (z_s[i] == depface[iface])
            vp[iface, 1] = alpha_s[i - 1]
            vs[iface, 1] = beta_s[i - 1]
            den[iface, 1] = rho[i - 1]
            vp[iface, 2] = alpha_s[i]
            vs[iface, 2] = beta_s[i]
            den[iface, 2] = rho[i]
            iflag[i - 1] = iface
            break
            println("ERROR: interface $iface not found")
        end
    end
end

for iwave in 1:2
    iw::Int64 = iwave
    pmin::Float64 = 1.25E-05 # avoid problems with flat earth near p=0
    
    if iwave == 1
        vel0::Float64 = vpsource
        pmax = 1. / (vpmin + 0.0001)
    else
        vel0 = vssource
        pmax = 1. / (vsmin + 0.0001)
    end

    dp::Float64 = (pmax - pmin) / float(np - 1)
    p1::Float64 = pmin
    p2::Float64 = pmin + dp
    vel0minus::Float64 = vel0 - 0.01 # used in takeoff angle calculation

    dang0::Float64 = asin(p2 * vel0minus) - asin(p1 * vel0minus)
    for i in 1:np
        amp0[i, iw] = 0
        p[i, iw] = pmin + float(i - 1) * dp
        if p[i, iw] * vel0 > 1.0 continue end
        if i == 1
            angcor::Float64 = 1.0
        else
            dang::Float64 = asin(p[i, iw] * vel0minus) - asin(p[i - 1, iw] * vel0minus)
            angcor = dang / dang0
        end
        sinthe = p[i, iw] * vel0minus
        amp0[i, iw] = sqrt(angcor * sinthe)
    end

    imth::Int64 = 3 # best for spherical earth
    for i in 1:npts-1
        h::Float64 = z[i+1] - z[i]
        ilay::Int64 = i
        for ip in 1:np
            dx[ilay, ip, iw], dt[ilay, ip, iw], irtr = LAYERTRACE(p[ip, iw], h, slow[i, iw],
                slow[i+1, iw], imth)
            if irtr == -1
                iddir[ilay, ip, iw] = 1
            elseif irtr == 1
                iddir[ilay, ip, iw] = 1
            elseif irtr == 0
                iddir[ilay, ip, iw] = -1
            elseif irtr == 2
                iddir[ilay, ip, iw] = -1
                dx[ilay, ip, iw] = 2. * dx[ilay, ip, iw]
                dt[ilay, ip, iw] = 2. * dt[ilay, ip, iw]
            else
                print("***ERROR: irtr = ", irtr, " ", i, " ", ilay, " ", ip, " ", p[ip, iw])
                break
            end

            dtstar[ilay, ip, iw] = dt[ilay, ip, iw] / q[i+1, iw]
            if iflag[i] != 0
                iface = iflag[i]
                r::Float64 = erad - depface[iface]
                flatfact::Float64 = r / erad
                pcor::Float64 = p[ip, iw] / flatfact
            # line 668
            end
        end
    end
end

# line 565