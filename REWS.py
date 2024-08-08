# This document contains a Python implementation of an analytical expression for the
# evaluation of the rotor-averaged deficit of a turbine placed in a Gaussian wake.
# The Gaussian wake is defined by:
#
#   - sigma:
#       The standard deviation of the Gaussian wake. If the wake source is a yawed
#       turbine, then sigma is the standard deviation of the vertical direction.
#
#   - xi:
#       The eccentricity of the elliptic contours of the wake which is related to
#       the standard deviations of the wake in the vertical (sigma_z) and lateral
#       (sigma_y) directions via
#
#           xi^2 = 1 - (sigma_y / sigma_z)^2
#
#   - omega:
#       is a wind-veer coefficient. The Gaussian wake that is subject to wind veer
#       is defined as 
#   
#           W = C * exp( -(y + omega * z)^2 / (2 * sigma_y^2)^2 ) * exp( -z^2 / (2 * sigma_z^2)^2 )
#       
#       The wind-veer coefficient omega is related to the difference in wind direction 
#       across the top and bottom tips of the turbine (da)
#
#           omega = da * x / D
#
#       where x is the streamwise distance between the considered turbine and the wake 
#       source, and D is the diameter of the wake source.
#
#   - rho, delta:
#       are the polar distances (radial and azimuthal) between the center of the wake 
#       and the center of the considered turbine. rho is normalised by the wake 
#       standard deviation sigma.
#
#   - gamma:
#       
#           gamma = radius / sigma
#
#       the ratio between the radius of the considered turbine and the standard
#       deviation of the wake sigma.

import numpy as np
from scipy.special import i0, i1
import math

def anl_rews(gamma, rho, delta, xi, omega):

    # The equation for the rotor-averaged deficit is spreadout across
    # many lines to avoid multiple evaluations of the same quantity.
    
    # gamma is the ratio between the radius of the considered turbine and the standard
    # deviation of the wake.
    g2 = gamma**2

    # Terms in the form "I_n(a/rho) / rho" are used, where I_n is the modified
    # Bessel function of order n. When rho = 0 (i.e., no lateral offset between
    # the wake center and the considered turbine), this limit is finite. 
    # To avoid any numerical issues, if rho = 0, we set rho to a small number
    if rho == 0: rho = 0.0001

    # Quantities related to xi and omega that are frequently used so we 
    # evaluate them once.
    xi2 = xi**2
    o2 = omega**2
    one_m_xi2 = 1 - xi2
    half_o2_m_xi2 = (o2 - xi2) / 2.0

    # Trignometric function of the angle delta
    cdelta = np.cos(delta)                  # cosine
    sdelta = np.sin(delta)                  # sine
    tdelta = sdelta / cdelta                # tan
    
    # New introduced angles 
    phi_ns = math.atan2(2.0 * omega, xi2 - o2)
    phi_s = math.atan2(omega + one_m_xi2 * tdelta / (omega * tdelta + 1), 1)
    phi = 2 * phi_s - phi_ns
    
    # New introduced length scales.
    # They are naturally normalised by the wake standard deviation sigma.
    # These length scales are squared (*_sq)
    sp_sq = one_m_xi2 / ( 1.0 + half_o2_m_xi2 )
    sns_sq = one_m_xi2 / (half_o2_m_xi2**2 + o2)**(0.5)
    ss_sq = one_m_xi2 * np.cos(phi_s) / (cdelta + omega * sdelta)

    # Some parameters used multiple times
    aux1 = sp_sq * rho
    aux2 = rho / ss_sq

    # Chi squared
    chi2 = aux1 * aux2 * sp_sq / 4.0 / sns_sq / ss_sq
    
    # The quantity P_ns
    # The trignometrics of the angle phi are used only once and hence are evaluated inline.
    pns = np.exp(-chi2 * np.cos(phi)) * np.cos(chi2 * np.sin(phi)) - 1.0
    
    # The quantity lambda / rho
    lam_rho = gamma * ss_sq / aux1
    
    # The argument of the Bessel functions
    kappa = gamma * aux2

    # The Bessel functions
    ii0 = i0(kappa)
    ii1 = i1(kappa)
    ### We can use recursion to get I_2 in terms of I_0 and I_1.
    ii2 = ii0 - 2.0 / kappa * ii1   
    
    # A quantity that used multiple times
    R_ss2 = 0.5 * g2 / sp_sq
    exp_R_ss2 = np.exp(-R_ss2)

    # The integral mu_0
    # Put the maximum number of iterations to be 20. Typical convergence is achieved
    # after 6-10 iterations.
    # The tolerance for convergence is set to 0.1%
    mu0 = sp_sq / g2 * exp_R_ss2 * Psi(R_ss2, ii0, kappa*ii1, aux2**2 * sp_sq, 20, 0.001)

    # Half the integral M_eta in the manuscript
    M = mu0 * (1.0 + 2.0 * pns) - np.exp(-R_ss2) / R_ss2 * pns * (lam_rho * ii1 + lam_rho * lam_rho * ii2)
    
    # The constant \hat{sigma}
    shat_sq = 1.0 / sp_sq + np.cos(2.0 * delta - phi_ns) / sns_sq
    
    # Return W/C: the normalised deficit referenced to the scaling function C(x)
    return 2 * M * np.exp(-rho * rho / 2.0 * shat_sq) 

def Psi(x, i00, ki11, tau_sq, ns, tol):
    # ns: is the maximum number of iterations to have a converged psi
    # tol: If this tolerance is achieved, iterations are stopped.
    fkm1 = 1
    gkm1 = 0
    sum0 = 0
    sum1 = 0
    for kk in range(ns):
        k = kk + 1
        fk = (fkm1 + tau_sq * gkm1) / k
        gk = (fk + 2 * gkm1) / (2 * k)
        rr = x**k
        new0 = fk * rr
        new1 = gk * rr
        sum0 += new0
        sum1 += new1
        if (abs(new0 / sum0) <= tol) and (abs(new1 / sum1) <= tol): break
        fkm1 = fk
        gkm1 = gk
    return i00 * sum0 - ki11 * sum1
