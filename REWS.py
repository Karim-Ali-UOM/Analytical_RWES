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
#
#   - n:
#       is the averaging order of the deficit by evaluating
#           W_n = [ (1 / A) * area_integral(W^n dA) ] ^ (1/n)
#       where A is the area of the disk representing the turbine

from scipy.special import i0, i1, owens_t, erf
from math import sqrt, sin, cos, pi, sqrt, atan2, atan, exp

def anl_rews(gamma, rho, delta, xi, omega, n = 1, method = "rectangular"):
    if method == "circular":
        rews = anl_rews_circular(gamma, rho, delta, xi, omega, n = n)
    elif method == "rectangular":
        rews = anl_rews_rectangular(gamma, rho, delta, xi, omega, n = n)
    return rews
    
def anl_rews_rectangular(gamma_in, beta_in, delta, xi, omega, ly = 0.88623, lz = 0.88623, n = 1):
    
    # ly: the ratio between the side length of the rectangular disk in y direction and the
    #     the actual radius of the turbine
    
    # lz: the ratio between the side length of the rectangular disk in z direction and the
    #     the actual radius of the turbine
    
    # n: averaging order
    
    # The ratio between the radius of the considered turbine and the standard
    # deviation of the wake.
    gamma = gamma_in * sqrt(n)
    
    # The ratio between the offset between the two turbines and the standard
    # deviation of the wake.
    rho = beta_in * sqrt(n)
    
    # Sine the angle delta
    sd = sin(delta)
    
    # Cosine the angle delta
    cd = cos(delta)
    
    # The ratio between the side length of the rectangular disk in y direction
    # and the standard deviation of the wake.
    gz = gamma * ly
    
    # The ratio between the side length of the rectangular disk in z direction
    # and the standard deviation of the wake.
    gy = gamma * lz
    
    # Save this value as it is used multiple times
    sqrt_om_xi2 = sqrt(1-xi*xi)
    
    # The ratio between the effect of wind veer and elliptic stretching due to
    # yaw misalignment of the wake source
    a = omega / sqrt_om_xi2
    
    # The generalised Owen T function is calculated at the four
    # vertices of the rectangular disk
    rhosd = rho * sd
    rhocd = rho * cd
    t = + gen_owen_t(rhosd - gz, a, (rhocd + gy) / sqrt_om_xi2) \
        - gen_owen_t(rhosd + gz, a, (rhocd + gy) / sqrt_om_xi2) \
        - gen_owen_t(rhosd - gz, a, (rhocd - gy) / sqrt_om_xi2) \
        + gen_owen_t(rhosd + gz, a, (rhocd - gy) / sqrt_om_xi2)
    return (pi * sqrt_om_xi2 / (2 * gy * gz) * t)**(1 / n)
    
def anl_rews_circular(gamma, rho, delta, xi, omega, n = 1):

    # The equation for the rotor-averaged deficit is spreadout across
    # many lines to avoid multiple evaluations of the same quantity.
    
    # n: averaging order
    
    # The ratio between the radius of the considered turbine and the standard
    # deviation of the wake.
    g2 = gamma**2

    # Terms in the form "I_n(a/rho) / rho" are used, where I_n is the modified
    # Bessel function of order n. When rho = 0 (i.e., no lateral offset between
    # the wake center and the considered turbine), this limit is finite. 
    # To avoid any numerical issues, if rho = 0, we set rho to a small number
    if rho == 0: rho = 0.0001
    
    # To avoid dividing by zero
    if xi == 0: xi = 0.0001
    if omega == 0: omega = 0.0001

    # Quantities related to xi and omega that are frequently used so we 
    # evaluate them once.
    xi2 = xi**2
    o2 = omega**2
    one_m_xi2 = 1 - xi2
    half_o2_m_xi2 = (o2 - xi2) / 2.0

    # Trignometric function of the angle delta
    cdelta = cos(delta)                  # cosine
    sdelta = sin(delta)                  # sine
    tdelta = sdelta / cdelta                # tan
    
    # New introduced angles 
    phi_ns = atan2(2.0 * omega, xi2 - o2)
    phi_s = atan2(omega + one_m_xi2 * tdelta / (omega * tdelta + 1), 1)
    phi = 2 * phi_s - phi_ns
    
    # New introduced length scales.
    # They are naturally normalised by the wake standard deviation sigma.
    # These length scales are squared (*_sq)
    sp_sq = one_m_xi2 / ( 1.0 + half_o2_m_xi2 ) / n
    sns_sq = one_m_xi2 / (half_o2_m_xi2**2 + o2)**(0.5) / n
    ss_sq = one_m_xi2 * cos(phi_s) / (cdelta + omega * sdelta) / n

    # Some parameters used multiple times
    aux1 = sp_sq * rho
    aux2 = rho / ss_sq

    # Chi squared
    chi2 = aux1 * aux2 * sp_sq / 4.0 / sns_sq / ss_sq
    
    # The quantity P_ns
    # The trignometrics of the angle phi are used only once and hence are evaluated inline.
    pns = exp(-chi2 * cos(phi)) * cos(chi2 * sin(phi)) - 1.0
    
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
    R2_sp2 = 0.5 * g2 / sp_sq
    exp_R2_sp2 = exp(-R2_sp2)

    # The integral mu_0
    # Put the maximum number of iterations to be 20. Typical convergence is achieved
    # after 6-10 iterations.
    # The tolerance for convergence is set to 0.1%
    mu0 = sp_sq / g2 * exp_R2_sp2 * Psi(R2_sp2, ii0, kappa*ii1, aux2**2 * sp_sq, 20, 0.001)
    
    # Half the integral M_eta in the manuscript
    M = mu0 * (1.0 + 2.0 * pns) - exp_R2_sp2 / R2_sp2 * pns * (lam_rho * ii1 + lam_rho * lam_rho * ii2)
    
    # The constant \hat{sigma}
    shat_sq = 1.0 / sp_sq + cos(2.0 * delta - phi_ns) / sns_sq
    
    # Return W/C: the normalised deficit referenced to the scaling function C(x)
    return (2 * M * exp(-rho * rho / 2.0 * shat_sq))**(1.0/n) 

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

def gen_owen_t(h,a,b):
    # Safe check against dividing by 0
    if(b==0): b = 0.001
    if(h==0): h = 0.001
        
    q1 = a + b/h
    q2 = (h + a*b + a*a*h) / b
    return -(atan(q1) + atan(q2)) / (2 * pi) + \
        owens_t(h, q1) + owens_t(b / sqrt(1 + a * a), q2)
