# --------------------------------------------------
# FORCE_RT forces RT coefficients to identical values
#   iforce=1 equalizes all scattered waves (0.25 for P/SV if all present)
#   iforce=2 for no conversions, equal scat. waves of same type
#
#  Note:  Any input rt that are zero are kept at zero
#         Remaining rt are scaled to still sum to 1
#

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