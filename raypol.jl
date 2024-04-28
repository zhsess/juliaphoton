function TO_CAR(the::Float64, phi::Float64, r::Float64)
    degrad::Float64 = pi / 180.0
    z::Float64 = r * cos(the * degrad)
    x::Float64 = r * sin(the * degrad) * cos(phi * degrad)
    y::Float64 = r * sin(the * degrad) * sin(phi * degrad)
    return x, y, z
end

function TO_POL(x::Float64, y::Float64, z::Float64)
    degrad::Float64 = pi / 180.0
    r::Float64 = sqrt(x ^ 2 + y ^ 2 + z ^ 2)
    the::Float64 = acos(z / r) * degrad
    if (x != 0.0) || (y != 0.0)
        phi::Float64 = atan(y, x) * degrad
    else
        phi = 0.0
    end
    return the, phi, r
end

function CROSS(v1::Array{Float64, 1}, v2::Array{Float64, 1})
    v3 = zeros(Float64, 3)
    v3[1] = v1[2] * v2[3] - v1[3] * v2[2]
    v3[2] = v1[3] * v2[1] - v1[1] * v2[3]
    v3[3] = v1[1] * v2[2] - v1[2] * v2[1]
    return v3
end

function VDOT(v1::Array{Float64, 1}, v2::Array{Float64, 1})
    dot::Float64 = v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]
    return dot
end

function ORTHOCHECK(v1::Array{Float64, 1}, v2::Array{Float64, 1}, v3::Array{Float64, 1})
    dot12 = VDOT(v1, v2)
    dot13 = VDOT(v1, v3)
    dot23 = VDOT(v2, v3)
    d1::Float64 = sqrt(v1[1] ^ 2 + v1[2] ^ 2 + v1[3] ^ 2)
    d2::Float64 = sqrt(v2[1] ^ 2 + v2[2] ^ 2 + v2[3] ^ 2)
    d3::Float64 = sqrt(v3[1] ^ 2 + v3[2] ^ 2 + v3[3] ^ 2)
    iret::Int64 = 0
    if abs(dot12) > 0.0001 iret = 1 end
    if abs(dot13) > 0.0001 iret = 1 end
    if abs(dot23) > 0.0001 iret = 1 end
    if abs(d1 - 1.0) > 0.0001 iret = 1 end
    if abs(d2 - 1.0) > 0.0001 iret = 1 end
    if abs(d3 - 1.0) > 0.0001 iret = 1 end
    return iret
end

# BENDRAY updates the ray and polarization direction vectors for
# the scattered ray.  The a1,a2,a3 vectors are overwritten with
# their new values.
# Inputs:
#     a1 = x1 axis in Sato and Fehler (S polarization direction)
#     a2 = x2 axis in S&F
#     a3 = x3 axis in S&F (ray direction)
#     psi = scattering angle from a3 (radians)
#     zeta = scattering angle from a1 (radians)
#     spol = scattered S polarization (radians) 
#            (=0 for pure psi polarization)
#            Note that this polarization is w.r.t the plane defined
#            by the incident and scattered ray vectors.  It is not
#            w.r.t the initial S polarization direction.
#
function BENDRAY(a1::Array{Float64, 1}, a2::Array{Float64, 1}, a3::Array{Float64, 1}, psi::Float64, zeta::Float64, spol::Float64)
    b1 = zeros(Float64, 3)
    b2 = zeros(Float64, 3)
    b3 = zeros(Float64, 3)
    c1 = zeros(Float64, 3)
    c2 = zeros(Float64, 3)
    c3 = zeros(Float64, 3)

    s_psi::Float64 = 0
    s_zeta::Float64 = 0
    i::Int64 = 0

    # define b vectors to be cartesian coor. of e_r, e_psi, e_zeta
    b3[3] = cos(psi)
    b3[1] = sin(psi) * cos(zeta)
    b3[2] = sin(psi) * sin(zeta)
    b1[3] = -sin(psi)
    b1[1] = cos(psi) * cos(zeta)
    b1[2] = cos(psi) * sin(zeta)
    b2 = CROSS(b3, b1)

    # define c vectors in terms of S polarization coordinate system
    s_psi = cos(spol)
    s_zeta = sin(spol)
    c1[1] = b1[1] * s_psi + b2[1] * s_zeta
    c1[2] = b1[2] * s_psi + b2[2] * s_zeta
    c1[3] = b1[3] * s_psi + b2[3] * s_zeta
    c3 = b3
    c2 = CROSS(c3, c1)

    # set b vectors to new ray based vectors using a and c vectors
    for i in 1:3
        b1[i] = c1[1] * a1[i] + c1[2] * a2[i] + c1[3] * a3[i]
        b2[i] = c2[1] * a1[i] + c2[2] * a2[i] + c2[3] * a3[i]
        b3[i] = c3[1] * a1[i] + c3[2] * a2[i] + c3[3] * a3[i]
    end

    a1 = b1
    a2 = b2
    a3 = b3

    return a1, a2, a3
end

# SCATRAYPOL computes change in ray parameter, ray azimuth,
# and vertical direction for 3-D scattering for both P, S 
# and P/S conversions.  This version includes S polarization.
#    Inputs:  p  =  table of ray parameters for both P and S
#            np  =  number of p values in p array
#            iw1 =  incident wave type (1=P, 2=S)
#            svsh = incident S polarization (radians, measured from SV pol.)
#            iw2 =  scattered wave type (1=P, 2=S)
#            ip  =  ray index for incident wave
#            psi =  scattering angle from ray direction (radians)
#            zeta = scattering angle from x1 axis (initial S polarization dir.)
#            spol = scattered S polarization (radians, 0=pure psi)
#            vp0  =  P velocity at scattering point
#            vs0  =  S velocity at scattering point
#            idir =  incident ray up/down direction (1=downgoing, 2=upgoing)
#            azi  =  azimuth of incident ray (degrees)
#   Returns: ip   =  ray index for scattered wave
#            idir =  scattered wave up/down direction
#            azi  =  azimuth of scattered wave (degrees)
#            svsh =  scattered S polarization (radians)
#
function SCATRAYPOL(p::Array{Float64, 2}, np::Int64, iw1::Int64, svsh::Float64, iw2::Int64, ip::Int64, psi::Float64, zeta::Float64, spol::Float64, vp0::Float64, vs0::Float64, idir::Int64, azi::Float64)
    nray0::Int64 = 50000

    vel0 = zeros(Float64, 2)
    slow = zeros(Float64, 2)
    a1 = zeros(Float64, 3)
    a2 = zeros(Float64, 3)
    a3 = zeros(Float64, 3)
    b1 = zeros(Float64, 3)
    b2 = zeros(Float64, 3)
    b3 = zeros(Float64, 3)

    degrad::Float64 = 180.0 / pi

    vel0[1] = vp0
    vel0[2] = vs0
    slow[1] = 1.0 / vp0
    slow[2] = 1.0 / vs0

    azi0::Float64 = azi
    svsh0::Float64 = svsh

    p1::Float64 = p[ip, iw1]
    if slow[iw1] < p1
        println("*** ERROR in SCATRAYP")
        idir = 0
        return ip, idir, azi, svsh
    end

    the::Float64 = asin(p1 * vel0[iw1]) * degrad # incident ray angle from vert (deg.)
    if idir == 1
        the = 180.0 - the
    end
    phi::Float64 = azi # incident angle from North

    a3[1], a3[2], a3[3] = TO_CAR(the, phi, 1.0)
    the2::Float64 = abs(the - 90.0)
    phi2::Float64 = phi
    if the < 90.0 phi2 = phi + 180.0 end

    a1[1], a1[2], a1[3] = TO_CAR(the2, phi2, 1.0)
    a2 = CROSS(a3, a1)

    iret = ORTHOCHECK(a1, a2, a3)
    if iret != 0
        println("*** ERROR in SCATRAYP: ORTHOCHECK failed")
    end

    for i in 1:3
        b3[i] = a3[i]
        b1[i] = cos(svsh) * a1[i] + sin(svsh) * a2[i]
        b2[i] = -sin(svsh) * a1[i] + cos(svsh) * a2[i]
    end

    b1, b2, b3 = BENDRAY(b1, b2, b3, psi, zeta, spol)
    the, phi, r = TO_POL(b3[1], b3[2], b3[3])

    p2::Float64 = slow[iw2] * sin(the / degrad)
    if the < 90.0
        idir = 2
    else
        idir = 1
    end
    azi = phi

    if round(azi/1000.0) != 0
        println("**SCATRAYPOL problem with azi ", azi)
        exit()
    end

    the2 = the - 90.0
    phi2 = phi
    if the < 90.0 phi2 = phi + 180.0 end
    a1[1], a1[2], a1[3] = TO_CAR(the2, phi2, 1.0)
    dot = VDOT(a1, b1)
    if dot > 1.0 dot = 1.0 end
    if dot < -1.0 dot = -1.0 end
    svsh = acos(dot)

    if round(svsh/10.0) != 0
        println("**SCATRAYPOL problem with svsh ", svsh)
        exit()
    end

    ptarg::Float64 = p2
    frac::Float64 = (ptarg - p[1, iw2]) / (p[np, iw2] - p[1, iw2])
    if frac < 0.0 frac = 0.0 end

    if (ip <= 0) || (ip > np + 20)
        println("**SCATRAYPOL problem with ip ", ip)
        exit()
    end

    if p[ip, iw2] > slow[iw2] ip = ip - 1 end
    if p[ip, iw2] > slow[iw2] ip = ip - 1 end
    if (ip < 1) || (ip > np)
        println("**WARNING")
    end

    return ip, idir, azi, svsh

end