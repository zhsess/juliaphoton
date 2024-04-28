# --------------------------------------------------
# FORCE_RT forces RT coefficients to identical values
#   iforce=1 equalizes all scattered waves (0.25 for P/SV if all present)
#   iforce=2 for no conversions, equal scat. waves of same type
#
#  Note:  Any input rt that are zero are kept at zero
#         Remaining rt are scaled to still sum to 1
#
# function FORCE_RT(rt::Array{Float64})
#     for i in 1:2 # P/Sv
#         for k in 1:2

#     end
# end
function FORCE_RT(rt::Array{Float64}, iforce::Int64)
    for i in 1:2
        for k in 1:2
            if iforce == 2
                rt[i, 1, k, 3-k] = 0.0
                rt[i, 2, k, 3-k] = 0.0
            end
            n::Int64 = 0
            for j in 1:2
                for l in 1:2
                    if rt[i, j, k, l] != 0 n = n + 1 end
                end
            end
            if n != 0
                rtfix::Float64 = 1.0 / float(n)
            else
                rtfix = 0.0
            end
            for j in 1:2
                for l in 1:2
                    if rt[i, j, k, l] != 0
                        rt[i, j, k, l] = rtfix
                    end
                end
            end
        end
    end

    for i in 1:2
        k::Int64 = 3
        n::Int64 = 0
        for j in 1:2
            l::Int64 = 3
            if rt[i, j, k, l] != 0 n = n + 1 end
        end
        if n != 0
            rtfix = 1.0 / float(n)
        else
            rtfix = 0.0
        end
        for j in 1:2
            l = 3
            if rt[i, j, k, l] != 0
                rt[i, j, k, l] = rtfix
            end
        end
    end
    return rt
end

# --------------------------------------------------
# subroutine RTCOEF calculates reflection/transmission coefficients
# for interface between two solid layers, based on the equations on 
# p. 150 of Aki and Richards.
#
#  Inputs:    vp1     =  P-wave velocity of layer 1 (top layer)
#  (real)     vs1     =  S-wave velocity of layer 1
#             den1    =  density of layer 1
#             vp2     =  P-wave velocity of layer 2 (bottom layer)
#             vs2     =  S-wave velocity of layer 2
#             den2    =  density of layer 2
#             hslow   =  horizontal slowness (ray parameter)
#  Returns:   rt(1)   =  down P to P up     (refl)
#  (complex)  rt(2)   =  down P to S up     (refl)
#             rt(3)   =  down P to P down   (tran)
#             rt(4)   =  down P to S down   (tran)
#             rt(5)   =  down S to P up     (refl)
#             rt(6)   =  down S to S up     (refl)
#             rt(7)   =  down S to P down   (tran)
#             rt(8)   =  down S to S down   (tran)
#             rt(9)   =    up P to P up     (tran)
#             rt(10)  =    up P to S up     (tran)
#             rt(11)  =    up P to P down   (refl)
#             rt(12)  =    up P to S down   (refl)
#             rt(13)  =    up S to P up     (tran)
#             rt(14)  =    up S to S up     (tran)
#             rt(15)  =    up S to P down   (refl)
#             rt(16)  =    up S to S down   (refl)
#
# NOTE:  All input variables are real.  
#        All output variables are complex!
#        Coefficients are not energy normalized.
#
function RTCOEF(vp1::Float64, vs1::Float64, den1::Float64, vp2::Float64, vs2::Float64, den2::Float64, hslow::Float64)

    alpha1::Complex{Float64} = complex(vp1, 0.0)
    beta1::Complex{Float64} = complex(vs1, 0.0)
    rho1::Complex{Float64} = complex(den1, 0.0)
    alpha2::Complex{Float64} = complex(vp2, 0.0)
    beta2::Complex{Float64} = complex(vs2, 0.0)
    rho2::Complex{Float64} = complex(den2, 0.0)
    p::Complex{Float64} = complex(hslow, 0.0)

    cone::Complex{Float64} = complex(1.0, 0.0)
    ctwo::Complex{Float64} = complex(2.0, 0.0)

    si1::Complex{Float64} = alpha1 * p
    si2::Complex{Float64} = alpha2 * p
    sj1::Complex{Float64} = beta1 * p
    sj2::Complex{Float64} = beta2 * p

    ci1::Complex{Float64} = sqrt(cone - si1 ^ 2)
    ci2::Complex{Float64} = sqrt(cone - si2 ^ 2)
    cj1::Complex{Float64} = sqrt(cone - sj1 ^ 2)
    cj2::Complex{Float64} = sqrt(cone - sj2 ^ 2)

    term1::Complex{Float64} = (cone - ctwo * beta2 * beta2 * p * p)
    term2::Complex{Float64} = (cone - ctwo * beta1 * beta1 * p * p)
    a::Complex{Float64} = rho2 * term1 - rho1 * term2
    b::Complex{Float64} = rho2 * term1 + ctwo * rho1 * beta1 * beta1 * p * p
    c::Complex{Float64} = rho1 * term2 + ctwo * rho2 * beta2 * beta2 * p * p
    d::Complex{Float64} = ctwo * (rho2 * beta2 * beta2 - rho1 * beta1 * beta1)
    e::Complex{Float64} = b * ci1 / alpha1 + c * ci2 / alpha2
    f::Complex{Float64} = b * cj1 / beta1 + c * cj2 / beta2
    g::Complex{Float64} = a - d * ci1 * cj2 / (alpha1 * beta2)
    h::Complex{Float64} = a - d * ci2 * cj1 / (alpha2 * beta1)
    den::Complex{Float64} = e * f + g * h * p * p

    rt = zeros(Complex{Float64}, 16)

    trm1::Complex{Float64} = b * ci1 / alpha1 - c * ci2 / alpha2
    trm2::Complex{Float64} = a + d * ci1 * cj2 / (alpha1 * beta2)
    rt[1] = (trm1 * f - trm2 * h * p * p) / den # refl down P to P up

    trm1 = a * b + c * d * ci2 * cj2 / (alpha2 * beta2)
    rt[2] = (- ctwo * ci1 * trm1 * p) / (beta1 * den) # refl down P to S up
    rt[3] = ctwo * rho1 * ci1 * f / (alpha2 * den) # trans down P to P down
    rt[4] = ctwo * rho1 * ci1 *h * p / (beta2 * den) # trans down P to S down

    trm1 = a * b + c * d * ci2 * cj2 / (alpha2 * beta2)
    rt[5] = (- ctwo * cj1 * trm1 * p) / (alpha1 * den) # refl down S to P up

    trm1 = b * cj1 / beta1- c * cj2 / beta2
    trm2 = a + d * ci2 * cj1 / (alpha2 * beta1)
    rt[6] = (trm1 * e - trm2 * g * p * p) / den # refl down S to S up
    rt[7] = - ctwo * rho1 * cj1 * g * p / (alpha2 * den) # trans down S to P down
    rt[8] = ctwo * rho1 * cj1 * e / (beta2 * den) # trans down S to S down

    trm1 = b * ci1 / alpha1 - c * ci2 / alpha2
    trm2 = a + d * ci2 * cj1 / (alpha2 * beta1)
    rt[11] = - (trm1 * f + trm2 * g * p * p) / den # refl up P to P down

    trm1 = a * c + b * d * ci1 / (alpha1 * beta1)
    rt[12] = (ctwo * ci2 * trm1 * p) / (beta2 * den) # refl up P to S down
    rt[9] = ctwo * rho2 * ci2 * f / (alpha1 * den) # trans up P to P up
    rt[10] = - ctwo * rho2 * ci2 * g * p / (beta1 * den) # trans up P to S up

    trm1 = a * c + b * d * ci1 * cj1 / (alpha1 * beta1)
    rt[15] = (ctwo * cj2 * trm1 * p) / (alpha2 * den) # refl up S to P down

    trm1 = b * cj1 / beta1 - c * cj2 / beta2
    trm2 = a + d * ci1 * cj2 / (alpha1 * beta2)

    rt[16] = (trm1 * e + trm2 * h * p * p) / den # refl up S to S down
    rt[13] = ctwo * rho2 * cj2 * h * p / (alpha1 * den) # trans up S to P up
    rt[14] = ctwo * rho2 * cj2 * e / (beta1 * den) # trans up S to S up

    return rt
end

# --------------------------------------------------
# subroutine RTCOEF_SH calculates SH reflection/transmission coefficients
# for interface between two solid layers, based on the equations on 
# p. 144 of Aki and Richards.
#
#  Inputs:    vs1     =  S-wave velocity of layer 1 (top layer)
#  (real)     den1    =  density of layer 1
#             vs2     =  S-wave velocity of layer 2
#             den2    =  density of layer 2
#             hslow   =  horizontal slowness (ray parameter)
#  Returns:   rt(1)   =  down S to S up     (refl)
#  (complex)  rt(2)   =  down S to S down   (tran)
#             rt(3)  =   up S to S up       (tran)
#             rt(4)  =   up S to S down     (refl)
#
# NOTE:  All input variables are real.  
#        All output variables are complex!
#        Coefficients are not energy normalized.
#
function RTCOEF_SH(vs1::Float64, den1::Float64, vs2::Float64, den2::Float64, hslow::Float64)

    rt = zeros(Complex{Float64}, 4)

    beta1::Complex{Float64} = complex(vs1, 0.0)
    rho1::Complex{Float64} = complex(den1, 0.0)
    beta2::Complex{Float64} = complex(vs2, 0.0)
    rho2::Complex{Float64} = complex(den2, 0.0)
    p::Complex{Float64} = complex(hslow, 0.0)

    cone::Complex{Float64} = complex(1.0, 0.0)

    sj1::Complex{Float64} = beta1 * p
    sj2::Complex{Float64} = beta2 * p
    cj1::Complex{Float64} = sqrt(cone - sj1 ^ 2)
    cj2::Complex{Float64} = sqrt(cone - sj2 ^ 2)

    D::Complex{Float64} = rho1 * beta1 * cj1 + rho2 * beta2 * cj2
    rt[1] = (rho1 * beta1 * cj1 - rho2 * beta2 * cj2) / D
    rt[4] = - rt[1]
    rt[2] = 2.0 * rho1 * beta1 * cj1 / D
    rt[3] = 2.0 * rho2 * beta2 * cj2 / D

    return rt
end

# --------------------------------------------------
# subroutine RTCOEF_POW calculates reflection/transmission coefficients
# for interface between two solid layers, based on the equations on 
# power as a real number.
#
#  Inputs:    vp1     =  P-wave velocity of layer 1 (top layer)
#  (real)     vs1     =  S-wave velocity of layer 1
#             den1    =  density of layer 1
#             vp2     =  P-wave velocity of layer 2 (bottom layer)
#             vs2     =  S-wave velocity of layer 2
#             den2    =  density of layer 2
#             p       =  horizontal slowness (ray parameter)
#  Returns:   refcoef =  2x2x3x3 array with coefficients
#                        i = incident wave direction (1=down, 2=up)
#                        j = scattered wave direction (1=down, 2=up)
#                        k = incident wave type (1=P, 2=SV, 3=SH)
#                        l = scattered wave type (1=P, 2=SV, 3=SH)
#
#  Requires:  RTCOEF, RTCOEF_SH
function RTCOEF_POW(vp1::Float64, vs1::Float64, den1::Float64, vp2::Float64, vs2::Float64, den2::Float64, p::Float64)
    v = zeros(Float64, 2, 3)
    den = zeros(Float64, 2)
    refcoef = zeros(Float64, 2, 2, 3, 3)

    v[1, 1] = vp1
    v[1, 2] = vs1
    v[1, 3] = vs1
    den[1] = den1
    v[2, 1] = vp2
    v[2, 2] = vs2
    v[2, 3] = vs2
    den[2] = den2

    rt = RTCOEF(vp1, vs1, den1, vp2, vs2, den2, p)
    rt_sh = RTCOEF_SH(vs1, den1, vs2, den2, p)

    refcoef[1, 2, 1, 1] = abs(rt[1])
    refcoef[1, 2, 1, 2] = abs(rt[2])
    refcoef[1, 1, 1, 1] = abs(rt[3])
    refcoef[1, 1, 1, 2] = abs(rt[4])
    refcoef[1, 2, 2, 1] = abs(rt[5])
    refcoef[1, 2, 2, 2] = abs(rt[6])
    refcoef[1, 1, 2, 1] = abs(rt[7])
    refcoef[1, 1, 2, 2] = abs(rt[8])
    refcoef[2, 2, 1, 1] = abs(rt[9])
    refcoef[2, 2, 1, 2] = abs(rt[10])
    refcoef[2, 1, 1, 1] = abs(rt[11])
    refcoef[2, 1, 1, 2] = abs(rt[12])
    refcoef[2, 2, 2, 1] = abs(rt[13])
    refcoef[2, 2, 2, 2] = abs(rt[14])
    refcoef[2, 1, 2, 1] = abs(rt[15])
    refcoef[2, 1, 2, 2] = abs(rt[16])

    refcoef[1, 2, 3, 3] = abs(rt_sh[1])
    refcoef[1, 1, 3, 3] = abs(rt_sh[2])
    refcoef[2, 2, 3, 3] = abs(rt_sh[3])
    refcoef[2, 1, 3, 3] = abs(rt_sh[4])

    for idir1 in 1:2
        ilay1::Int64 = idir1
        for idir2 in 1:2
            ilay2::Int64 = 3 - idir2
            for iw1 in 1:3
                for iw2 in 1:3
                    if refcoef[idir1, idir2, iw1, iw2] == 0.0
                        continue
                    end
                    slow1::Float64 = 1.0 / v[ilay1, iw1]
                    slow2::Float64 = 1.0 / v[ilay2, iw2]
                    if (p >= slow1) || (p >= slow2)
                        refcoef[idir1, idir2, iw1, iw2] = 0.0
                        continue
                    end

                    cos1::Float64 = sqrt(1.0 - p ^ 2 * v[ilay1, iw1] ^ 2)
                    cos2::Float64 = sqrt(1.0 - p ^ 2 * v[ilay2, iw2] ^ 2)
                    f1::Float64 = v[ilay1, iw1] * den[ilay1] * cos1
                    f2::Float64 = v[ilay2, iw2] * den[ilay2] * cos2

                    rtnorm::Float64 = refcoef[idir1, idir2, iw1, iw2]
                    rtnorm = (f2 / f1) * rtnorm ^ 2
                    refcoef[idir1, idir2, iw1, iw2] = sqrt(rtnorm)
                end
            end
        end
    end

    # the following section is a kluge to eliminate transmitted
    # energy at surfaces for the case where we have used very small
    # velocities to simulate air or fluid.

    for idir1 in 1:2
        ilay1 = idir1
        for idir2 in 1:2
            ilay2 = 3 - idir2
            if v[ilay2, 1] < 0.2
                refcoef[idir1, idir2, 1, 1] = 0.0
                refcoef[idir1, idir2, 2, 1] = 0.0
            end
            if v[ilay2, 2] < 0.1
                refcoef[idir1, idir2, 1, 2] = 0.0
                refcoef[idir1, idir2, 2, 2] = 0.0
            end
            if v[ilay2, 3] < 0.1
                refcoef[idir1, idir2, 3, 3] = 0.0
            end
        end
    end
    # now renormalize P/SV so sum of energy is one
    for i in 1:2
        for k in 1:2
            sum::Float64 = 0.0
            for j in 1:2
                for l in 1:2
                    sum = sum + refcoef[i, j, k, l]
                end
            end
            if sum == 0.0
                continue
            end
            for j in 1:2
                for l in 1:2
                    refcoef[i, j, k, l] = refcoef[i, j, k, l] / sum
                end
            end
        end
    end
    # now renormalize SH so sum of energy is one
    for i in 1:2
        k::Int64 = 3
        sum = 0.0
        for j in 1:2
            l::Int64 = 3
            sum = sum + refcoef[i, j, k, l]
        end
        if sum == 0.0
            continue
        end
        for j in 1:2
            l = 3
            refcoef[i, j, k, l] = refcoef[i, j, k, l] / sum
        end
    end
    return refcoef
end