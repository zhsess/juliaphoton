# calculate flat earth transformation
function FLATTEN(z_s::Float64, vel_s::Float64)
    erad:: Float64 = 6371.0
    r::Float64 = erad - z_s
    z_f::Float64 = - erad * log(r / erad)
    vel_f::Float64 = vel_s * (erad / r)
    return z_f, vel_f
end