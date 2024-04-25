# LAYERTRACE calculates the travel time and range offset
# for ray tracing through a single layer.
#
# Input:    p     =  horizontal slowness
#           h     =  layer thickness
#           utop  =  slowness at top of layer
#           ubot  =  slowness at bottom of layer
#           imth  =  interpolation method
#                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)
#                         = 2,  v(z) = a - b*z
#                         = 3,  v(z) = a*exp(-b*z)
#
# Returns:  dx    =  range offset
#           dt    =  travel time
#           irtr  =  return code
#                 = -1, zero thickness layer
#                 =  0,  ray turned above layer
#                 =  1,  ray passed through layer
#                 =  2,  ray turned within layer, 1 segment counted
#
# Note:  This version does calculation in double precision,
#        but all i/o is still single precision
#
function LAYERTRACE(p::Float64, h::Float64, utop::Float64, ubot::Float64, imth::Int64)
    dx::Float64 = 0.0
    dt::Float64 = 0.0
    irtr::Int64 = 0

    if h == 0.0 # check for zero thickness layer
        dx = 0.0
        dt = 0.0
        irtr = -1
        return dx, dt, irtr
    end

    u::Float64 = utop
    y::Float64 = u - p
    if y <= 0.0 # complex vertical slowness
        dx = 0.0
        dt = 0.0
        irtr = 0
        return dx, dt, irtr
    end

    q::Float64 = y * (u + p)
    qs::Float64 = sqrt(q)
    # special function needed for integral at top of layer
    if imth == 2
        y = u + qs
        if p != 0.0 y = y / p end
        qr = log(y)
    elseif imth == 3
        qr = atan(qs, p)
    end

    if imth == 1
        b::Float64 = - (utop ^ 2 - ubot ^ 2) / (2.0 * h)
    elseif imth == 2
        vtop::Float64 = 1.0 / utop
        vbot::Float64 = 1.0 / ubot
        b = (vbot - vtop) / h
    else
        b = log(utop / ubot) / h
    end

    if b == 0.0 # constant slowness layer
        b = 1.0 / h
        etau::Float64 = qs
        ex::Float64 = p / qs
        irtr = 1
        dx = ex / b
        dtau::Float64 = etau / b
        dt = dtau + p * dx
        return dx, dt, irtr
    end

    # integral at upper limit, 1/b factor omitted until end
    if imth == 1
        etau = - q * qs / 3.0
        ex = -qs * p
    elseif imth == 2
        ex = qs / u
        etau = qr - ex
        if p != 0.0 ex = ex / p end
    else
        etau = qs - p * qr
        ex = qr
    end

    # check lower limit to see if we have turning point
    u = ubot
    if u <= p # if turning point, no contribution from lower layers
        irtr = 0
        dx = ex / b
        dtau = etau / b
        dt = dtau + p * dx
        return dx, dt, irtr
    end

    irtr = 1
    q = (u - p) * (u + p)
    qs = sqrt(q)
    if imth == 1
        etau = etau + q * qs / 3.0
        ex = ex + qs * p
    elseif imth == 2
        y = u + qs
        z::Float64 = qs / u
        etau = etau = z
        if p != 0.0
            y = y / p
            z = z / p
        end
        qr = log(y)
        etau = etau - qr
        ex = ex - z
    else
        qr = atan(qs, p)
        etau = etau - qs + p * qr
        ex = ex - qr
    end
    dx = ex / b
    dtau = etau / b
    dt = dtau + p * dx
    return dx, dt, irtr
end