# Program to model seismic phases as energy packets
# Originally written 2003-2004 by Peter Shearer (pshearer@ucsd.edu)
# Converted to f90 in 2012 by Nick Mancinelli (n.j.mancinelli@gmail.com)
# Converted to Julia in 2024 by Hao Zhang (zhsess@gmail.com)
# Last modified: 2024-04-24

using SpecialFunctions
include("gavg.jl")
include("layertrace.jl")
include("rt.jl")
include("utilities.jl")
include("raypol.jl")

# read in parameters
vmodel::String = "./iasp91.smo"
zsource::Float64 = 0.0 # assumed source depth (km)
freq::Float64 = 1.0 # assumed freq. (Hz), affects attenuation and scattering
om::Float64 = 2.0 * pi * freq

iforcerefl::Int64 = 0 # 0=normal, 1=force equal refl/trans

xwind1::Float64 = 0.0 # xwind1,xwind2,twind1,twind2 for ray info dump
xwind2::Float64 = 0.0
twind1::Float64 = 0.0
twind2::Float64 = 0.0

iwstart::Int64 = 7 # radiate: (1) P,  (2) SH, (3) SV, (4) SH/SV, (5) P+S, (6) custom Es/Ep, (7) custom pmax
pmax_takeoff::Float64 = 0.07 # maximum takeoff p (custom pmax option 7 above only)

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
nscatmin::Int64 = 0        # min,max number of scattering events for output
nscatmax::Int64 = 999999

idebug::Int64 = 2 # out.debug: 0=none, 1=full, 2=scat
if idebug != 0
    debugfile::String = "OUTPUT/out.debug"
    f_debug = open(debugfile, "w")
    f_surface = open("OUTPUT/out.surface", "w")
end
idebug0::Int64 = idebug
nraydump::Int64 = 50000 # write output at multiples of this number of rays

# if iforcerefl = -1
nslowout::Int64 = 5
nslowdim::Int64 = 49

slowstack = zeros(Float64, nslowdim, nslowdim, nslowout)
xslow1 = zeros(Float64, nslowout)
xslow2 = zeros(Float64, nslowout)
tslow1 = zeros(Float64, nslowout)
tslow2 = zeros(Float64, nslowout)

# define variables
nface0::Int64 = 6
nlay0::Int64 = 655
nray0::Int64 = 50000
erad::Float64 = 6371.0
ecircum::Float64 = 2.0 * pi * erad
kmdeg::Float64 = ecircum / 360.0
degrad::Float64 = 180.0 / pi

ntdim::Int64 = 3001
nxdim::Int64 = 360

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

rt = zeros(Float64, 2, 2, 3, 3, nface0, nray0, 2)
depsave = zeros(Float64, 10000, 2)
iwsave = zeros(Int64, 1000)

mutable struct SavedValues
    firstcall::BitArray
    psi2::Array{Float64, 2}
    zeta2::Array{Float64, 2}
    sumg::Array{Float64, 3}
    nk::Int64
end
save_data= SavedValues([true, true, true, true, true, true], 
    zeros(Float64, 66000, 6),
    zeros(Float64, 66000, 6),
    zeros(Float64, 66000, 4, 6),
    0)

rbin = zeros(Float64, ntdim, nxdim)
rbin_p = zeros(Float64, ntdim, nxdim)
rbin_sv = zeros(Float64, ntdim, nxdim)
rbin_sh = zeros(Float64, ntdim, nxdim)
rbin_z = zeros(Float64, ntdim, nxdim)
rbin_rad = zeros(Float64, ntdim, nxdim)
rbin_zcore = zeros(Float64, ntdim, nxdim)

stnmin::Float64 = 0.0
stnmax::Float64 = 0.0
dummy = zeros(Float64, 10)
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

    println("Scattering Strengths for layer $iscat:")
    println("gpp0 = $gpp0, gps0 = $gps0, gsp0 = $gsp0, gss0 = $gss0")
    println("afp = ", 1. / (gpp0 + gps0), ", ", 1. / (gsp0 + gss0))
    println("gppm = $gppm, gpsm = $gpsm, gspm = $gspm, gssm = $gssm")
    println("moment afp =", 1. / (gppm + gpsm), ", ", 1. / (gspm + gssm), "\n")

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
        global ilay::Int64 = i
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

                rt[1:2, 1:2, 1:3, 1:3, iface, ip, iw] = RTCOEF_POW(vp[iface, 1], vs[iface, 1], den[iface, 1], 
                    vp[iface, 2], vs[iface, 2], den[iface, 2], pcor)
                
                if iforcerefl == 1
                    rt[1:2, 1:2, 1:3, 1:3, iface, ip, iw] = FORCE_RT(rt[1:2, 1:2, 1:3, 1:3, iface, ip, iw], iforcerefl)
                end

                if (idebug == 3) && (iw == 1) && (iface == 2)
                    println("$iw, $ip, $(p[ip, iw]), $(rt[1,1,1,1,iface,ip,iw]), $(rt[2,1,1,1,iface,ip,iw])")
                end
            end
        end
    end
    nlay::Int64 = ilay
    println("Finished ray tracing, (iwave, nlay, np) = ($iw, $nlay, $np) \n")
end

# ----------- start spraying rays ------------
tmax::Float64 = 3000.0
nray::Int64 = 0
nsurf::Int64 = 0
if (iwstart == 1) || (iwstart == 7)
    iww1::Int64 = 1
    iww2::Int64 = 1
elseif iwstart <= 4
    iww1 = 2
    iww2 = 2
else
    iww1 = 1
    iww2 = 2
end

if izsource == 1
    idirbeg::Int64 = 1
    idirend::Int64 = 1
else
    idirbeg = -1
    idirend = 1
end

nsweep::Int64 = 0
while nsweep < 100 # loop over all rays (150)
    println(nsweep)
    for iw_init in iww1:iww2 # loop over wave types (800)
        for idir_init in idirbeg:2:idirend # loop over directions (780)
            global nsweep += 1
            for ipp in 1:np # loop over ray parameters (760)
                ifinish::Bool = true

                ip::Int64 = ipp
                iw::Int64 = iw_init
                if (iwstart == 7) && (p[ip, iw] > pmax_takeoff)
                    continue
                end

                idir::Int64 = idir_init
                ampstart::Float64 = amp0[ip, iw]
                if ampstart == 0.0 continue end
                svsh::Float64 = 0.0

                if iw == 2
                    if iwstart != 6
                        ampstart = ampstart * sqrt(23.4)
                    else
                        ampstart = ampstart * sqrt(spenergy_ratio)
                    end

                    svsh = 0.0
                    if iwstart == 2
                        svsh = pi / 2.0
                    elseif (iwstart >= 4) && (iwstart < 7)
                        fran::Float64 = rand()
                        svsh = fran * pi
                    end
                end

                global nray = nray + 1
                if nray % 100 == 0 # record debug for every 100 rays
                    global idebug = idebug0
                else
                    global idebug = 0
                end

                if (idebug == 1) || (idebug == 3)
                    write(f_debug, "SPRAY: $nray $ipp $(p[ip,iw]) $iw $idir \n")
                end

                ixdir::Int64 = 1
                if idir == 1
                    global ilay = izsource
                else
                    global ilay = izsource - 1
                end

                # initialize ray variables
                x::Float64 = 0.0
                t::Float64 = 0.0
                plat::Float64 = 0.0
                plon::Float64 = 0.0
                azi::Float64 = 90.0
                tstar::Float64 = 0.0
                isave::Int64 = 0
                nscat::Int64 = 0
                iscatold::Int64 = -99

                while true # 200
                    x = x + dx[ilay, ip, iw] * float(ixdir) # in km
                    t = t + dt[ilay, ip, iw] * float(ixdir) # in seconds
                    if t > tmax
                        break
                    end
                    
                    tstar = tstar + dtstar[ilay, ip, iw]
                    idir = idir * iddir[ilay, ip, iw]

                    if idebug == 1
                        write(f_debug, "$nray $ip $iw $idir $ilay $t $(z_s[ilay]) \n")
                    end

                    # check for volume scattering
                    if (iflagvol[ilay] != 0) && (dt[ilay, ip, iw] != 0)
                        iscat::Int64 = iflagvol[ilay]
                        # to allow stronger near-source scat.
                        if (xmaxscat[iscat] < 180.0)
                            xdeg::Float64 = x / kmdeg
                            plat2, plon2 = SPH_LOC(plat, plon, xdeg, azi)
                            xkm::Float64 = (90.0 - plat2) * kmdeg
                            if (x > xmaxscat[iscat]) iscat = iscat - 1 end
                            if iscat < 1
                                println("**iscat sync problem: &iscat, $ilay, $x \n")
                                break
                            end
                        end
    
                        if iscat != iscatold # entered new scattering volume
                            global afp::Float64 = 1.0 / scatprob[iscat, iw, 1]
                            global freepath::Float64 = -log(rand()) * afp
                            global distsum::Float64 = 0.0
                        end
                        iscatold = iscat
    
                        if iw == 1
                            dist::Float64 = dt[ilay, ip, iw] * alpha[ilay]
                        else
                            dist = dt[ilay, ip, iw] * beta[ilay]
                        end
                        distsum = distsum + dist

                        # ----- scattering section -----
                        if distsum >= freepath # scattered
                            nscat = nscat + 1
                            iscatold = -99 # reset 

                            if idebug != 0
                                write(f_debug, "SC: $nray $ip $iw $(p[ip, iw]) $idir $ilay $iscat \n")
                                write(f_debug, "    $nray $nscat $distsum $freepath $afp $x $t \n")
                            end

                            distsum = 0.0
                            iwave1::Int64 = iw
                            iwoff::Int64 = (iw - 1) * 2 # for RANDANG2 input

                            fran = rand() * scatprob[iscat, iw, 1]

                            if idebug != 0
                                write(f_debug, "    $fran $(scatprob[iscat, iw, 2]) \n")
                            end

                            if fran <= scatprob[iscat, iw, 2]
                                psi2, zeta2, spol = RANDANG2(iscat, iwoff+1, el[iscat], nu[iscat], gam0[iscat], eps[iscat], alen[iscat], save_data)
                                iwave2::Int64 = 1
                            else
                                psi2, zeta2, spol = RANDANG2(iscat, iwoff+2, el[iscat], nu[iscat], gam0[iscat], eps[iscat], alen[iscat], save_data)
                                iwave2 = 2
                            end
    
                            volang::Float64 = psi2 * degrad # input angle for subfunctions
                            
                            if idebug == 1
                                write(f_debug, "Scat: $iw, $iwave2, $volang \n")
                            end
                            if idebug != 0
                                write(f_debug, "A: $volang, $psi2, $zeta2, $spol, $iwave2 \n")
                            end
    
                            xdeg = x / kmdeg
                            plat2, plon2 = SPH_LOC(plat, plon, xdeg, azi)
    
                            # following codes fix numerical problems near podes
                            if xdeg < 0.5
                                azi2 = azi
                            elseif (xdeg > 179.5) && (xdeg < 180.5)
                                azi2 = 180.0 - azi
                                if azi2 < 0.0 azi2 = azi2 + 360.0 end
                            else
                                del, azi2 = SPH_AZI(plat2, plon2, plat, plon)
                                azi2 = azi2 + 180.0
                                if azi2 > 360.0 azi2 = azi2 - 360.0 end
                            end
    
                            if idebug == 2
                                write(f_debug, "L: $plat, $plon, $azi, $xdeg \n")
                                write(f_debug, "   $plat2, $plon2, $azi2 \n")
                            end

                            # azimuth calculation has problems near poles
                            # the following check should exit to next ray if NaN occurs
    
                            if isnan(azi2) then
                                println("ERROR: azimuth is nan, skip to next ray")
                                println("azi2 = $azi2")
                                ifinish = false
                                break
                            end
    
                            plat = plat2
                            plon = plon2
                            azi = azi2
                            x = 0.0
                            distsum = 0.0
                            
                            # scat point at end of ray path
                            idir2::Int64 = 0
                            while idir2 == 0
                                if idir == 1
                                    vp0::Float64 = alpha[ilay+1]
                                    vs0::Float64 = beta[ilay+1]
                                else
                                    vp0 = alpha[ilay]
                                    vs0 = beta[ilay]
                                end
                                idir2 = idir
                                azi2 = azi
                                
                                if idebug != 0
                                    write(f_debug, "InSCATRAYPOL: $np, $iwave1, $svsh, $iwave2, $ip \n")
                                    write(f_debug, "$psi2, $zeta2, $spol, $vp0, $vs0 \n")
                                    write(f_debug, "$idir2, $azi2 \n")
                                end

                                ip, idir2, azi2, svsh = SCATRAYPOL(p, np, iwave1, svsh, iwave2, ip, psi2, zeta2, spol, vp0, vs0, idir2, azi2)
                                
                                if idebug != 0
                                    write(f_debug, "OutSCATRAYPOL: $np, $iwave1, $svsh, $iwave2, $ip \n")
                                    write(f_debug, "$psi2, $zeta2, $spol, $vp0, $vs0 \n")
                                    write(f_debug, "$idir2, $azi2 \n")
                                    write(f_debug, "S: $iw, $ipp, $ip, $azi, $azi2, $plat, $plon, $idir \n")
                                    write(f_debug, "\n")
                                end

                                if idir2 == 0
                                    println("ERROR: idir2 = 0")
                                    idir = -1
                                end
                            end

                            azi = azi2
                            iwave = iwave2
                            iw = iwave
                            if idir2 != idir
                                if idir2 == -1 & iddir[ilay, ip, iw] == -1
                                    idir2 = 1
                                end
                                idir = idir2
                                continue
                            end
                        end
                    end

                    # Interface section
                    if idir == 1 # downgoing
                        if iflag[ilay + 1] == 0
                            ilay = ilay + 1
                        else
                            iface = iflag[ilay + 1]
                            # for dumping ray info
                            isave = isave + 1
                            depsave[isave, 1] = depface[iface]
                            iwsave[isave] = iw
                            
                            if iw == 1 # incident P wave
                                fran = rand()
                                rt1::Float64 = rt[1, 1, iw, iw, iface, ip, iw] # down trans. no conversion
                                rt2::Float64 = rt[1, 2, iw, iw, iface, ip, iw] # reflected, no conversion
                                rt3::Float64 = rt[1, 1, iw, 3 - iw, iface, ip, iw] # down trans, phase conver.
                                rt4::Float64 = rt[1, 2, iw, 3 - iw, iface, ip, iw] # reflected, phase conver.
                                test::Float64 = rt1 + rt2 + rt3 + rt4
                                if abs(test - 1.0) > 0.01
                                    println("***Warning in downgoing ref/trans coef ")
                                    println("(iw, iface, ip, depface[iface]) = ($iw, $iface, $ip, $(depface[iface]))")
                                    println("rt1, rt2, rt3, rt4 = $rt1, $rt2, $rt3, $rt4")
                                    r::Float64 = erad - depface[iface]
                                    flatfact::Float64 = r / erad
                                    pcor = p[ip, iw] / flatfact
                                    println("$(p[ip, iw]), $pcor, $(1.0/(vp[iface, 1])), $(1.0/(vp[iface, 2]))")
                                    xdeg = x / kmdeg
                                    tmin::Float64 = t / 60.
                                    println("xdeg, tmin = $xdeg, $tmin")
                                    rt2 = 1.0
                                end
                                test1::Float64 = rt1
                                test2::Float64 = rt2 + test1
                                test3::Float64 = rt3 + test2

                                if idebug == 1
                                    write(f_debug, "Pdown face: $ip, $iw, $(depface[iface]) \n")
                                end

                                if fran <= test1
                                    ilay = ilay + 2
                                elseif fran <= test2
                                    idir = -1
                                elseif fran <= test3
                                    ilay = ilay + 2
                                    iw2 = 3 - iw
                                    ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                    if ip2 == 0
                                        println("***DEBUG2: $ilay, $fran, $test1, $test2, $test3")
                                        exit()
                                    end
                                    iw = iw2
                                    ip = ip2
                                    svsh = 0.0
                                else
                                    idir = -1
                                    iw2 = 3 - iw
                                    ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                    if ip2 == 0
                                        println("***DEBUG3: $ilay, $fran, $test1, $test2, $test3")
                                        exit()
                                    end
                                    iw = iw2
                                    ip = ip2
                                    svsh = 0.0
                                end
                            else # incident S wave
                                svfrac::Float64 = cos(svsh) ^ 2
                                shfrac::Float64 = 1.0 - svfrac
                                test1 = svfrac * rt[1, 1, 2, 1, iface, ip, iw] # SV to P trans
                                test2 = test1 + svfrac * rt[1, 2, 2, 1, iface, ip, iw] # SV to P refl
                                test3 = test2 + svfrac * rt[1, 1, 2, 2, iface, ip, iw] + shfrac * rt[1, 1, 3, 3, iface, ip, iw] # SV to SH / SH to SV
                                fran = rand()

                                if idebug == 1
                                    write(f_debug, "Sdown face: $ip, $iw, $(depface[iface]) \n")
                                end

                                if fran <= test1 # SV to P trans.
                                    ilay = ilay + 2
                                    iw2 = 3 - iw
                                    ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                    if ip2 == 0
                                        println("***DEBUG3: $ilay, $fran, $test1, $test2, $test3")
                                        exit()
                                    end
                                    iw = iw2
                                    ip = ip2
                                    svsh = 0.0
                                elseif fran <= test2 # SV to P refl.
                                    idir = -1
                                    iw2 = 3 - iw
                                    ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                    if ip2 == 0
                                        println("***DEBUG4: $ilay, $fran, $test1, $test2, $test3")
                                        exit()
                                    end
                                    iw = iw2
                                    ip = ip2
                                    svsh = 0.0
                                elseif fran <= test3 # SV to SH / SH to SV trans.
                                    ilay = ilay + 2
                                    svamp::Float64 = sqrt(svfrac * rt[1, 1, 2, 2, iface, ip, iw])
                                    shamp::Float64 = sqrt(shfrac * rt[1, 1, 3, 3, iface, ip, iw])
                                    svsh = atan(shamp, svamp)
                                end
                            end
                        end
                    else # upgoing (idir = -1)
                        if ilay == 1 # at surface, need to output t, x
                            isave = isave + 1
                            depsave[isave, 1] = 0.0
                            depsave[isave, 2] = 0.0
                            iwsave[isave] = iw

                            tmin = t / 60.0
                            xdeg = x / kmdeg
                            iwrap::Int64 = mod(Int64(floor(xdeg / 180.0)), 2)
                            plat2, plon2 = SPH_LOC(plat, plon, xdeg, azi)
                            xdeg, azidum = SPH_AZI(0.0, 0.0, plat2, plon2)

                            psecdeg::Float64 = p[ip, iw] * 111.19
                            xdum, azidum1 = SPH_AZI(plat2, plon2, 0.0, 0.0)
                            xdum, azidum2 = SPH_AZI(plat2, plon2, plat, plon)
                            
                            if iwrap == 1 azidum2 = azidum2 + 180.0 end
                            slowang::Float64 = azidum2 - azidum1
                            psecdegtran::Float64 = sin(slowang / degrad) * psecdeg
                            psecdegrad::Float64 = cos(slowang / degrad) * psecdeg

                            if idebug == 1
                                write(f_debug, "at surface: $ip, $iw, $xdeg, $tmin \n")
                                write(f_surface, "$nray, $ip, $iw, $xdeg, $tmin \n")
                            end

                            ix::Int64 = round(xdeg * 2 + 0.5)
                            if ix < 1 ix = 1 end
                            if ix > nxdim ix = nxdim end
                            it::Int64 = round(t + 0.5)
                            if it < 1 it = 1 end
                            if it > ntdim it = ntdim end

                            amp::Float64 = 25.0 * ampstart * exp(-freq * pi * tstar) # attenuation
                            if amp < 1.0E-12 amp = 0.0 end # avoind underflow errors
                            amp2::Float64 = amp ^ 2

                            if (nscat >= nscatmin) && (nscat <= nscatmax)
                                rbin[it, ix] = rbin[it, ix] + amp2
                                if iw == 1
                                    rbin_p[it, ix] = rbin_p[it, ix] + amp2
                                    sinthe::Float64 = p[ip, iw] * vpmin
                                    amp_rad::Float64 = sinthe * sqrt(amp2)
                                    energy_rad::Float64 = amp_rad ^ 2
                                    energy_vert::Float64 = amp2 - energy_rad
                                    rbin_z[it, ix] = rbin_z[it, ix] + energy_vert
                                    if psecdeg <= 5.0
                                        rbin_zcore[it, ix] = rbin_zcore[it, ix] + energy_vert
                                    end
                                    rbin_rad[it, ix] = rbin_rad[it, ix] + energy_rad

                                    if iforcerefl == -1
                                        for k in 1:nslowout
                                            if (xdeg >= xslow1[k]) && (xdeg <= xslow2[k]) && (t >= tslow1[k]) && (t <= tslow2[k])
                                                i = (nslowdim + 1) / 2 + round(psecdegrad * 2.0)
                                                j = (nslowdim + 1) / 2 + round(psecdegtran * 2.0)
                                                if (i >= 1) && (i <= nslowdim) && (j >= 1) && (j <= nslowdim)
                                                    slowstack[i, j, k] = slowstack[i, j, k] + energy_vert
                                                end
                                            end
                                        end
                                    end
                                else
                                    svfrac = cos(svsh) ^ 2
                                    shfrac = 1.0 - svfrac
                                    rbin_sv[it, ix] = rbin_sv[it, ix] + amp2 * svfrac
                                    rbin_sh[it, ix] = rbin_sh[it, ix] + amp2 * shfrac
                                    sinthe = p[ip, iw] * vsmin
                                    amp_vert::Float64 = sinthe * sqrt(amp2)
                                    energy_vert = amp_vert ^ 2
                                    energy_rad = amp2 - energy_vert
                                    rbin_z[it, ix] = rbin_z[it, ix] + energy_vert
                                    if psecdeg <= 5.0
                                        rbin_zcore[it, ix] = rbin_zcore[it, ix] + energy_vert
                                    end
                                    rbin_rad[it, ix] = rbin_rad[it, ix] + energy_rad
                                end
                            end

                            if (xdeg >= xwind1) && (xdeg <= xwind2) && (t >= twind1) && (t <= twind2)
                                println("Match to window, xdeg, tmin, p:")
                                println("$xdeg, $tmin, $(p[ip, iw]), $tstar, $amp, $nscat")
                                nsave::Int64 = isave
                                for isave in 1:nsave
                                    println("$(depsave[isave, 1]), $(depsave[isave, 2]), $(iwsave[isave])")
                                end
                            end

                            global nsurf = nsurf + 1
                            if iw == 1 # upgoing P-wave at surface
                                fran = rand()
                                iface = 1
                                if iforcerefl != 2
                                    test1 = rt[2, 1, iw, iw, iface, ip, iw] # reflected no conversion
                                else
                                    test1 = 1.1
                                end

                                if fran <= test1
                                    idir = 1
                                else # phase conversion
                                    idir = 1
                                    iw2 = 3 - iw
                                
                                    ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                    if ip2 == 0
                                        println("***DEBUG3: $ilay, $fran, $test1")
                                        exit()
                                    end
                                    iw = iw2
                                    ip = ip2
                                    svsh = 0.0
                                end
                            else # upgoing S-wave at surface
                                svfrac = cos(svsh) ^ 2
                                shfrac = 1.0 - svfrac
                                iface = 1
                                test1 = svfrac * rt[2, 1, 2, 1, iface, ip, iw] # SV to P refl.
                                fran = rand()
                                if fran <= test1 # SV to P refl.
                                    idir = 1
                                    iw2 = 3 - iw
                                    ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                    if ip2 == 0
                                        println("***DEBUG3: $ilay, $fran, $test1, $test2, $test3")
                                        exit()
                                    end
                                    iw = iw2
                                    ip = ip2
                                    svsh = 0.0
                                else
                                    idir = 1 # S to S reflected
                                    svamp = sqrt(svfrac * rt[2, 1, 2, 2, iface, ip, iw]) # SV to SV refl.
                                    shamp = sqrt(shfrac * rt[2, 1, 3, 3, iface, ip, iw]) # SH to SH refl.
                                    svsh = atan(shamp, svamp)
                                end
                            end
                        else # upgoing, not at surface
                            if (ilay == 2) || (iflag[ilay - 1] == 0) # no interface
                                ilay = ilay - 1
                            else
                                iface = iflag[ilay - 1]
                                
                                isave = isave + 1 # for dumping ray info
                                depsave[isave, 1] = depface[iface]
                                depsave[isave, 2] = 0.0
                                iwsave[isave] = iw

                                if iw == 1 # incident P wave
                                    fran = rand()
                                    test1 = rt[2, 2, iw, iw, iface, ip, iw] # up trans. no conversion
                                    test2 = test1 + rt[2, 1, iw, iw, iface, ip, iw] # reflected, no conversion
                                    test3 = test2 + rt[2, 2, iw, 3 - iw, iface, ip, iw] # up trans., phase conver.

                                    if idebug == 1
                                        write(f_debug, "Pup face: $ip, $iw, $(depface[iface]) \n")
                                    elseif idebug == 3
                                        write(f_debug, "Pup face: $ip, $iw, $(depface[iface]) \n")
                                        write(f_debug, "    $(p[ip, iw]), $(rt[2, 1, iw, iw, iface, ip, iw]) \n")
                                    end

                                    if fran <= test1 # transmitted
                                        ilay = ilay - 2
                                    elseif fran <= test2 # reflected
                                        idir = 1
                                    elseif fran <= test3 # transmitted, phase conversion
                                        ilay = ilay - 2
                                        iw2 = 3 - iw
                                        ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                        if ip2 == 0
                                            println("***DEBUG6: $ilay, $fran, $test1, $test2, $test3")
                                            exit()
                                        end
                                        iw = iw2
                                        ip = ip2
                                        svsh = 0.0
                                    else # reflected, phase conversion
                                        idir = 1
                                        iw2 = 3 - iw
                                        ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                        if ip2 == 0
                                            println("***DEBUG7: $ilay, $fran, $test1, $test2, $test3")
                                            exit()
                                        end
                                        iw = iw2
                                        ip = ip2
                                        svsh = 0.0
                                    end
                                else # incident S wave
                                    svfrac = cos(svsh) ^ 2
                                    shfrac = 1.0 - svfrac
                                    test1 = svfrac * rt[2, 2, 2, 1, iface, ip, iw] # SV to P
                                    test2 = test1 + svfrac * rt[2, 1, 2, 1, iface, ip, iw] # SV to P refl.
                                    test3 = test2 + svfrac * rt[2, 2, 2, 2, iface, ip, iw] + shfrac * rt[2, 2, 3, 3, iface, ip, iw] # SV to SH / SH to SV

                                    fran = rand()

                                    if idebug == 1
                                        write(f_debug, "Sup face: $ip, $iw, $(depface[iface]) \n")
                                    end

                                    if fran <= test1 # SV to P trans
                                        ilay = ilay - 2
                                        iw2 = 3 - iw
                                        ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                        if ip2 == 0
                                            println("***DEBUG8: $ilay, $fran, $test1, $test2, $test3")
                                            exit()
                                        end
                                        iw = iw2
                                        ip = ip2
                                        svsh = 0.0
                                    elseif fran <= test2 # SV to P refl
                                        idir = 1
                                        iw2 = 3 - iw
                                        ip2 = UPDATE_P(p, np, iw, iw2, ip)
                                        if ip2 == 0
                                            println("***DEBUG9: $ilay, $fran, $test1, $test2, $test3")
                                            exit()
                                        end
                                        iw = iw2
                                        ip = ip2
                                        svsh = 0.0
                                    elseif fran <= test3 # SV to SH / SH to SV
                                        ilay = ilay - 2
                                        svamp = sqrt(svfrac * rt[2, 2, 2, 2, iface, ip, iw])
                                        shamp = sqrt(shfrac * rt[2, 2, 3, 3, iface, ip, iw])
                                        svsh = atan(shamp, svamp)
                                    else # S to S reflected
                                        idir = 1
                                        svamp = sqrt(svfrac * rt[2, 1, 2, 2, iface, ip, iw]) # SV to SV refl.
                                        shamp = sqrt(shfrac * rt[2, 1, 3, 3, iface, ip, iw]) # SH to SH refl.
                                        svsh = atan(shamp, svamp)
                                    end
                                end
                            end
                        end
                    end
                end

                if ifinish == false
                    continue
                end

                # if finished, write data
                if mod(nray, nraydump) == 0
                    println("Writing rbin files... nray, nsurf = $nray, $nsurf")
                    out_file::String = "out.photon"
                    f_out = open(out_file, "w+")
                    write(f_out, "0.0, 180.0, $nxdim, 0.0, 3000.0, $ntdim, $iwstart, $stnmin, $stnmax, $zsource, $zsource, $nray, $nsurf \n")
                    write(f_out, "$dummy \n")
                    for ix in 1:nxdim
                        for it in 1:ntdim
                            write(f_out, "$(rbin[it, ix]) ")
                        end
                    end
                    close(f_out)

                    out_file = "OUTPUT/out.photon_p"
                    f_out = open(out_file, "w+")
                    write(f_out, "0.0, 180.0, $nxdim, 0.0, 3000.0, $ntdim, $iwstart, $stnmin, $stnmax, $zsource, $zsource, $nray, $nsurf \n")
                    write(f_out, "$dummy \n")
                    for ix in 1:nxdim
                        for it in 1:ntdim
                            write(f_out, "$(rbin_p[it, ix]) ")
                        end
                    end
                    close(f_out)

                    out_file = "OUTPUT/out.photon_sv"
                    f_out = open(out_file, "w+")
                    write(f_out, "0.0, 180.0, $nxdim, 0.0, 3000.0, $ntdim, $iwstart, $stnmin, $stnmax, $zsource, $zsource, $nray, $nsurf \n")
                    write(f_out, "$dummy \n")
                    for ix in 1:nxdim
                        for it in 1:ntdim
                            write(f_out, "$(rbin_sv[it, ix]) ")
                        end
                    end
                    close(f_out)

                    out_file = "OUTPUT/out.photon_z"
                    f_out = open(out_file, "w+")
                    write(f_out, "0.0, 180.0, $nxdim, 0.0, 3000.0, $ntdim, $iwstart, $stnmin, $stnmax, $zsource, $zsource, $nray, $nsurf \n")
                    write(f_out, "$dummy \n")
                    for ix in 1:nxdim
                        for it in 1:ntdim
                            write(f_out, "$(rbin_z[it, ix]) ")
                        end
                    end
                    close(f_out)

                    out_file = "OUTPUT/out.photon_zcore"
                    f_out = open(out_file, "w+")
                    write(f_out, "0.0, 180.0, $nxdim, 0.0, 3000.0, $ntdim, $iwstart, $stnmin, $stnmax, $zsource, $zsource, $nray, $nsurf \n")
                    write(f_out, "$dummy \n")
                    for ix in 1:nxdim
                        for it in 1:ntdim
                            write(f_out, "$(rbin_zcore[it, ix]) ")
                        end
                    end
                    close(f_out)

                    out_file = "OUTPUT/out.photon_rad"
                    f_out = open(out_file, "w+")
                    write(f_out, "0.0, 180.0, $nxdim, 0.0, 3000.0, $ntdim, $iwstart, $stnmin, $stnmax, $zsource, $zsource, $nray, $nsurf \n")
                    write(f_out, "$dummy \n")
                    for ix in 1:nxdim
                        for it in 1:ntdim
                            write(f_out, "$(rbin_rad[it, ix]) ")
                        end
                    end
                    close(f_out)

                    println("Finished writing.")

                    n = 0
                    rsum::Float64 = 0.0
                    rmax::Float64 = 0.0
                    for it in 1:ntdim
                        for ix in 1:nxdim
                            if rbin[it, ix] == 0.0
                                continue
                            end
                            if rbin[it, ix] > rmax
                                rmax = rbin[it, ix]
                            end
                            n = n + 1
                            rsum = rsum + rbin[it, ix]
                        end
                    end
                    println("n, rmax, rsum = $n, $rmax, $rsum")
                end
            end
        end
    end
end

if idebug0 != 0
    close(f_debug)
    close(f_surface)
end
# line 565