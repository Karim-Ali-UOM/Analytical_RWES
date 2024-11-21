from numpy import linspace, pi, sqrt, meshgrid, zeros_like,	exp, mean
from numpy import array, zeros, zeros_like, power, divide, multiply
from numpy import cross, dot, clip, arctan2
from math import cos, sin, atan2, tan, atan, log
import matplotlib.patches as patches
from matplotlib.legend_handler import HandlerTuple
from scipy.special import i0, i1, iv, owens_t
from os.path import join
import time
import random
from scipy.interpolate import interp1d
from numpy.linalg import norm

class HandlerTupleVertical(HandlerTuple):
    """Plots all the given Lines vertical stacked."""

    def __init__(self, **kwargs):
        """Run Base Handler."""
        HandlerTuple.__init__(self, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        """Create artists (the symbol) for legend entry."""
        # How many lines are there.
        numlines = len(orig_handle)
        handler_map = legend.get_legend_handler_map()

        # divide the vertical space where the lines will go
        # into equal parts based on the number of lines
        height_y = (height / numlines)

        leglines = []
        for i, handle in enumerate(orig_handle):
            handler = legend.get_legend_handler(handler_map, handle)

            legline = handler.create_artists(legend, handle,
                                             xdescent,
                                             (2*i + 1)*height_y,
                                             width,
                                             2*height,
                                             fontsize, trans)
            leglines.extend(legline)
        return leglines

def wake_map():
    dia = 1.25
    thetas = linspace(0,2*pi,100,endpoint=True)
    wake_center = 0 + 0j
    turb_center = 1 + 1j
    sigma = 1.5
    xi = 0.4
    sigma_z = sigma
    sigma_y = sigma * sqrt(1-xi**2)
    chi = -0.6
    ys = linspace(-1.5,2,50)
    zs = linspace(-1.5,2,50)
    y, z = meshgrid(ys, zs, indexing='ij')
    wake = zeros_like(y)
    for i in range(y.shape[0]):
        for k in range(y.shape[1]):
            dy = y[i,k] - wake_center.real
            dz = z[i,k] - wake_center.imag
            wake[i,k] = exp(-(dy + chi*dz)**2/(2*sigma_y**2)) * exp(-dz**2/(2*sigma_z**2))
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif"
    })
    #fig = plt.figure(figsize=(3,3))
    #ax = plt.axes()
    fig, axs = plt.subplots(ncols = 2, nrows = 1, figsize=(6,3))
    ax = axs[0]
    ax2 = axs[1]
    
    im = ax2.contourf(y,z,wake, 300, cmap = "viridis", vmin=0, vmax=1)
    ax2.contour(y,z,wake, 17, colors = "grey", alpha = 0.7, linewidths = 0.5)
    for c in im.collections:
        c.set_edgecolor("face")
    #cb_ax = fig.add_axes([0.905, 0.0575, 0.02, 0.855])
    cb_ax = fig.add_axes([0.47, 0.06, 0.41, 0.02])
    cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal', format='%.1f')
    cbar.ax.tick_params(labelsize=6)
    #cbar.set_ticks(linspace(0,1,10))
    cbar.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    cbar.ax.set_title("$W / C$", fontsize = 6, x=1.05, y=1.0, pad=0)
    #AngleAnnotation((wake_center.real,wake_center.imag), (wake_center.real + 1,wake_center.imag), (0.8*wake_center.real+0.2*turb_center.real,0.8*wake_center.imag+0.2*turb_center.imag), ax=ax, size=75, unit="pixels", text="", linewidth = 0.5)
    
    for axx in axs:
        axx.spines['left'].set_visible(False)
        axx.spines['bottom'].set_visible(False)
        axx.spines['right'].set_visible(False)
        axx.spines['top'].set_visible(False)
        axx.set_yticks([])
        axx.set_xticks([])
        axx.set_aspect("equal")
    
    # upwind turbine
    r = dia / 2.0
    deflection = 0.2*dia
    upTurbine = [(r*cos(theta) + 1j * r * sin(theta)) + wake_center - deflection for theta in thetas]
    yup = [p.real for p in upTurbine]
    zup = [p.imag for p in upTurbine]
    for axx in axs:
        axx.plot(yup, zup, linewidth = 0.75, linestyle = "--", color = "tab:red")
    ax.scatter(wake_center.real - deflection, wake_center.imag, color = "tab:red", s = 5, zorder = 10)
    angg = 135*pi/180.0
    ax.plot([wake_center.real - deflection, wake_center.real - deflection + r * cos(angg)], [wake_center.real, wake_center.real + r * sin(angg)], linestyle = "-", color = "k", linewidth = 0.5, zorder = 1)
    ax.text(wake_center.real - deflection + r/2.0*cos(angg) - 0.15, wake_center.imag + r/2.0*sin(angg) - 0.1, "$R_o$", size = 7)
    downnTurbine = [(r*cos(theta) + 1j * r * sin(theta)) + turb_center for theta in thetas]
    yd = [p.real for p in downnTurbine]
    zd = [p.imag for p in downnTurbine]
    #ax.plot(yd, zd, linewidth = 0.75, linestyle = "-", color = "k", zorder = 10)
    for axx in axs:
        circle = plt.Circle((turb_center.real, turb_center.imag), dia/2, color='none', zorder = 100, ec = "tab:red")
        axx.add_patch(circle)
    ax.plot([wake_center.real, turb_center.real-wake_center.real], [wake_center.imag, wake_center.imag], linestyle = "-", color = "k", linewidth = 0.5, zorder = 1)
    ax.plot([turb_center.real, turb_center.real], [wake_center.imag, turb_center.imag], linestyle = "-", color = "k", linewidth = 0.5, zorder = 1)
    ax.plot([wake_center.real, turb_center.real], [wake_center.imag, turb_center.imag], linestyle = "-", color = "k", linewidth = 0.5, zorder = 1)
    ax.plot([wake_center.real - deflection, wake_center.real], [wake_center.imag, wake_center.imag], linestyle = "-", color = "k", linewidth = 0.5, zorder = 1)
    ang = 135 * pi / 180.0
    ax.plot([turb_center.real, turb_center.real + dia/2*cos(ang)], [turb_center.imag, turb_center.imag+dia/2*sin(ang)], linestyle = "-", color = "k", linewidth = 0.5, zorder = 1)
    ax.text(turb_center.real + dia/4*cos(ang) - 0.15, turb_center.imag + dia/4*sin(ang) - 0.1, "$R$", size = 7)
    for axx in axs:
       axx.text(wake_center.real + 0.22, wake_center.imag-0.15, "$y'$", size = 7)
       axx.text(wake_center.real - 0.05, wake_center.imag + 0.22 + 0.15, "$z'$", size = 7)
       axx.arrow(wake_center.real, wake_center.imag, 0.25, 0, color = "k", width = 0.01, zorder = 10)
       axx.arrow(wake_center.real, wake_center.imag, 0, 0.25, color = "k", width = 0.01, zorder = 10)
       axx.scatter(wake_center.real, wake_center.imag, color = "k", s = 10, zorder = 10)
    
    for axx in axs:
        axx.text(turb_center.real + 0.22, turb_center.imag-0.15, "$y$", size = 7)
        axx.text(turb_center.real - 0.05, turb_center.imag + 0.22 + 0.15, "$z$", size = 7)
        axx.scatter(turb_center.real, turb_center.imag, color = "k", s = 10, zorder = 10)
        axx.arrow(turb_center.real, turb_center.imag, 0.25, 0, color = "k", width = 0.01, zorder = 10)
        axx.arrow(turb_center.real, turb_center.imag, 0, 0.25, color = "k", width = 0.01, zorder = 10)
    
    ax.text(wake_center.real - deflection + 0.075, wake_center.imag + 0.05, "$d_o$", size = 7)
    ax.text((wake_center.real+turb_center.real)/2.0 + 0.1, wake_center.imag + 0.05, "$\Delta_y$", size = 7)
    ax.text(turb_center.real - 0.16, (wake_center.imag+turb_center.imag)/2 , "$\Delta_z$", size = 7)
    ax.text((0.6*turb_center.real + 0.4*wake_center.real), (0.6*wake_center.imag+0.4*turb_center.imag) +0.3, "$\\rho$", size = 7)
    
    radius = 0.16
    angd = 10*pi/180 + atan2(turb_center.imag-wake_center.imag, turb_center.real-wake_center.real)
    angd2 = -7.5*pi/180.0
    ts = linspace(0,angd, 10, endpoint=True)
    pxs = [radius*cos(t) for t in ts]
    pys = [radius*sin(t) for t in ts]
    # ax.plot(pxs, pys, linewidth = 0.5, color = "k")
    ax.text(wake_center.real+0.2, wake_center.imag+0.05, "$\delta$", size = 7)
    style = "Simple, tail_width=0.5, head_width=2.25, head_length=2.25"
    kw = dict(arrowstyle=style, color="k")
    a3 = patches.FancyArrowPatch((radius*cos(angd2), radius*sin(angd2)), (radius*cos(angd), radius*sin(angd)),connectionstyle="arc3,rad=0.4",
        linewidth=0.25, **kw)
    ax.add_patch(a3)
    
    ly = r*0.9
    lz = r*0.9
    rect = matplotlib.patches.Rectangle(xy=(turb_center.real-ly, turb_center.imag-lz),
            width = 2*ly, height = 2*lz, linestyle= "--", fc = "none", edgecolor="tab:blue")
    axs[0].add_patch(rect)
    rect2 = matplotlib.patches.Rectangle(xy=(turb_center.real-ly, turb_center.imag-lz),
            width = 2*ly, height = 2*lz, linestyle= "--", fc = "none", edgecolor="tab:blue")
    axs[1].add_patch(rect2)
    
    arrow_width = 0.05
    widtt = 0.0075 #0.01
    zrect = turb_center.imag + lz + 0.075 + 0.05
    yrect = turb_center.real - 0.075
    axs[0].text(yrect, zrect, "$2L_y$", size = 7)
    axs[0].arrow(yrect + 0.22, zrect + 0.02, 0.63*ly, 0,
        color = "k", width = widtt, zorder = 10,
        linewidth = arrow_width)
    axs[0].arrow(yrect - 0.05, zrect + 0.02, - 0.67*ly, 0,
        color = "k", width = widtt, zorder = 10,
        linewidth = arrow_width)
    
    
    axs[0].text(yrect+0.65*ly, zrect-0.75*lz, "$2L_z$", size = 7)
    axs[0].arrow(yrect+0.65*ly + 0.1, zrect-0.75*lz + 0.1 ,0, 0.27*lz,
        color = "k", width = widtt, zorder = 10,
        linewidth = arrow_width)
    axs[0].arrow(yrect+0.65*ly + 0.1, zrect-0.75*lz - 0.025 ,0, -1.35*lz,
        color = "k", width = widtt, zorder = 10,
        linewidth = arrow_width)
    
    axs[0].text(-1,1.95,"(a)",fontsize=8)
    axs[1].text(-1.5,2.1,"(b)",fontsize=8)
    fig.subplots_adjust(bottom=0.075, top=0.95, left=0.025, right=0.88, wspace=0.05, hspace=0.1)
    fig.savefig("wake_map.png", dpi=300)
    
def validate_one_case():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif"
    })
    plt.rc('text.latex', preamble='\\usepackage{braket}')

    fig, axs = plt.subplots(ncols = 4, nrows = 1, figsize=(7,3))
    daoa = 7*pi/180
    xds = [4,6,8,10]
    yaw = 20*pi/180 
    betas = [0.001,0.5,1,1.5,2,2.5,3,3.5,4]
    betas_refined = linspace(0.0001, 4, 50, endpoint=True)
    deltas = [0, pi/4, 3*pi/4]
    markers = ["o","^","s","v","D"]
    delta_handles = []
    delta_labels = []
    anl_handles = []
    num_handles = []
    colors = ["k", "tab:red", "tab:blue", "tab:orange"]
    letters = ["a","b","c","d","e","f","g","h","i","j",
        "k","l","m","n","o","p","q","r","s","t","u",
        "v","w","x","y","z"]
    delta_names = ["0", "\\pi/4", "3\\pi/4"]
    radius = 124.6 / 2.0 # NREL-5MW
    
    pss = fetch_points(radius, "sunflower", n = 2000)
    
    column_index = -1

    ct = 0.8
    ti = 0.05

    for xd in xds:
        gamma, xi, omega = calc_constants(ct, ti, yaw, xd, daoa)
        column_index += 1
        ax = axs[column_index]

        for delta in deltas:
            idelta = deltas.index(delta)
            label = "$\\delta=" + delta_names[idelta] + "$"
            
            curve,  = ax.plot(betas_refined, [solve_new(gamma, beta, xi, omega, delta) for beta in betas_refined],
                linewidth = 1, color = colors[idelta],
                label= label, zorder = 1)
            
            #curve,  = ax.plot(betas_refined, [solve_new_square(gamma, beta, xi, omega, delta) for beta in betas_refined],
            #    linewidth = 1, color = colors[idelta],
            #    label= label, zorder = 1)
            
            points = ax.scatter(betas, [solve_numerical(pss, radius, gamma, beta, delta, xi, omega) for beta in betas],
                marker=markers[idelta],
                s = 16,facecolors='white', edgecolor = colors[idelta],
                zorder = 20, label= label, linewidth = 0.75)

            #ax.plot(betas_refined, [solve_numerical([0], radius, gamma, beta, delta, xi, omega) for beta in betas_refined],
            #    linewidth = 0.75, color = colors[idelta],
            #    label= label, zorder = 1, linestyle="--")
            
            if column_index == 0:
                delta_handles.append((curve,points))
                delta_labels.append(label)

                anl_handles.append(curve)
                num_handles.append(points)
            ax.text(0.37,0.9,"(" + letters[column_index] + ") $x/D_o=" + str(xd) + "$",transform=ax.transAxes,fontsize = 10)
            ax.text(0.54,0.8,"$\\xi=" + str(round(xi,2)) + "$",transform=ax.transAxes,fontsize = 8)
            ax.text(0.54,0.73,"$\\omega=" + str(round(omega,2)) + "$",transform=ax.transAxes,fontsize = 8)
            ax.text(0.54,0.66,"$R/\\sigma=" + str(round(gamma,1)) + "$",transform=ax.transAxes,fontsize = 8)
            if idelta==0: ax.text(0.54,0.59,"$\\kappa^{(1)}=" + str(round(calc_kappa(gamma, xi, omega, delta),2)) + "$",
                c = "k",transform=ax.transAxes,fontsize = 8)
            #ax.text(0.05,0.25 - 0.07*idelta,"$\\mu_0^{(1)}=" + str(round(calc_mu0(gamma, xi, omega, delta, 1),2)) + "$",
            #    c = colors[idelta],transform=ax.transAxes,fontsize = 8)
            
            #ax.text(0.02,0.18,"$C_t=" + str(ct) + "$",transform=ax.transAxes,fontsize = 9)
            #ax.text(0.02,0.06,"$T_i=" + str(int(ti*100)) + "\%$",transform=ax.transAxes,fontsize = 9)
            
            ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
            if column_index == 0:
                ax.set_ylabel("$\\overline{W}_r^{(1)}/C$", fontsize = 10)
                #ax.set_ylabel("$\\braket{\\overline{W}^{(1)},\\hat{W}}/C$", fontsize = 9)
                ax.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
            else:
                ax.yaxis.set_ticklabels([])
            
            ax.set_xticks([0,1,2,3,4])
            ax.set_xlabel("$\\rho/\\sigma$", fontsize = 10)
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax.set_xlim([-0.2,4.2])
            ax.set_ylim([-0.03,0.7])   

    fig.legend(handles = delta_handles, labels = delta_labels, 
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.06, 0.92], ncol=3,
        handlelength = 2.5,
        prop={'size': 10}, frameon=False)
    fig.legend(handles = [(anl_handles[0],anl_handles[1],anl_handles[2])], labels = ["analytical"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.68, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    fig.legend(handles = [(num_handles[0],num_handles[1],num_handles[2])], labels = ["numerical"], 
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.83, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    
    fig.subplots_adjust(bottom=0.135, top=0.92, left=0.075, right=0.975, wspace=0.1, hspace=0.1)
    fig.savefig("validate_circle_single_7v_20y.png", dpi=500)

def solve_new_square(gammain, betain, xi, omega, delta, params = None, n=1):
    gamma = gammain* sqrt(n)
    beta = betain * sqrt(n)
    sd = sin(delta)
    cd = cos(delta)
    gz = gamma * 0.9
    gy = gamma * 0.9
    a = omega/sqrt(1-xi*xi)
    t = + gen_owen_t(beta*sd - gz, a, (beta*cd + gy)/sqrt(1-xi*xi)) \
        - gen_owen_t(beta*sd + gz, a, (beta*cd + gy)/sqrt(1-xi*xi)) \
        - gen_owen_t(beta*sd - gz, a, (beta*cd - gy)/sqrt(1-xi*xi)) \
        + gen_owen_t(beta*sd + gz, a, (beta*cd - gy)/sqrt(1-xi*xi))

    if t<0: t=0
    return (pi*sqrt(1-xi*xi) / (2*gy*gz) * t)**(1/n)

def gen_owen_t(h,a,b):
    q1 = a+b/h
    q2 = (h+a*b+a*a*h)/b
    qq = b/sqrt(1+a*a)
    g = -(atan(q1) + atan(q2))/(2*pi)
    return g + owens_t(h,q1) + owens_t(qq, h/b + a + a*a*h/b)

def calc_kappa(gamma, xi, omega, delta, n=1):
    sp_sq, sm_sq, sns_sq, phi_ns, ss_sq, phi_s = calc_states_4(xi, omega, delta)
    return n*gamma**2/(6*sns_sq)

def solve_new(gammain, betain, xi, omega, delta, n=1):
    gamma = gammain
    beta = betain

    sp_sq, _, sns_sq, phi_ns, ss_sq, phi_s = calc_states_4(xi, omega, delta)
    sp_sq = sp_sq / n
    sns_sq = sns_sq / n
    ss_sq = ss_sq / n
    
    chi = beta*sp_sq/(2*sqrt(sns_sq)*ss_sq)
    phi = 2*phi_s - phi_ns
    pns = exp(-chi**2 * cos(phi)) * cos(chi**2 * sin(phi)) - 1
    # sum of bessels
    summ = 0
    for nu in range(2):
        val = iv(nu+1,gamma*beta/ss_sq) * (gamma*ss_sq/(beta*sp_sq))**(nu+1)
        summ += val
    G = -2*sp_sq/gamma**2 * exp(-gamma**2/(2*sp_sq)) * pns * summ
    mu0 = sp_sq/gamma**2 * exp(-gamma**2/(2*sp_sq)) * generic_Phi(gamma, beta, 1.0/sqrt(ss_sq), 1.0/sqrt(sp_sq),100)
    G += mu0*(1+2*pns)
    if G < 0: G=0
    
    shat_sq = 1/(1/sp_sq + cos(2*delta-phi_ns)/sns_sq)
    return (2*exp(-beta**2/(2*shat_sq))*G)**(1/n)

def get_bc(alpha, n):
    ps, qs = get_pq(alpha, n)
    if alpha == 0:
        bs = [p for p in ps[1:] if abs(p) != 0]
        cs = zeros_like(bs)
    else:
        bs = [ps[2*j+2] * alpha for j in range(-1+int(n/2))]
        cs = [qs[2*j+3] for j in range(-1+int(n/2))]
    return bs, cs

def get_pq(alpha, n):
    p = [-1/alpha, 0]
    q = [0, 0]
    for k in range(n):
        if k in [0,1]: continue
        p.append((2 * p[k-2] - alpha * q[k-1]) / k)
        q.append((2 * q[k-2] + alpha * p[k-1])/(k-1))
    return p, q

def solve_numerical(pss, R, gamma, beta, delta, xi, omega, n=1):
    sigma = R / gamma
    sigma_z = sigma / sqrt(n)
    sigma_y = sigma * sqrt(1-xi**2) / sqrt(n)
    turb_center = beta*sigma * cos(delta) + 1j* beta*sigma * sin(delta)
    ps = [p + turb_center for p in pss]
    ws = [exp(-(p.real + omega*p.imag)**2 / 2.0 / sigma_y**2) * exp(-(p.imag)**2 / 2.0 / (sigma_z**2)) for p in ps]
    return (mean(ws))**(1/n)

def calc_states_4(xi, omega, delta):
    sp_sq = (1-xi**2) / (1+omega**2 / 2 - xi**2 / 2)
    sm_sq = (1-xi**2) / (-omega**2 / 2 + xi**2 / 2)
    sns_sq = (1-xi**2) * ((-omega**2 / 2 + xi**2 / 2)**2 + omega**2)**(-0.5)
    phi_ns = atan2(2*omega, xi**2-omega**2)
    ss_sq = (((cos(delta)+omega*sin(delta))/(1-xi**2))**2 + ((omega*cos(delta)+omega**2 *sin(delta))/(1-xi**2)+sin(delta))**2)**(-0.5)
    phi_s = atan2(omega+(1-xi**2)*tan(delta)/(omega*tan(delta) + 1), 1)
    return sp_sq, sm_sq, sns_sq, phi_ns, ss_sq, phi_s

def calc_constants(ct, ti, yaw, xd, daoa):
    eps = 0.2*sqrt((1+sqrt(1-ct))/(2*sqrt(1-ct)))
    ks = 0.003678 + 0.3837 * ti
    sz0 = sqrt((1+sqrt(1-ct*cos(yaw)))/(8*(1+sqrt(1-ct))))
    gamma = 0.5/(ks*xd+eps)
    omega = daoa * xd
    syd = ks*xd + cos(yaw)*sz0
    szd = ks*xd + sz0
    xi = sqrt(1-(syd/szd)**2)
    return gamma, xi, omega

def sunflower(n, alpha):
    b = round(alpha * sqrt(n))  # number of boundary points
    phi = (sqrt(5) + 1) / 2  # golden ratio
    rs= []
    ths = []
    for k in range(1, n + 1):
        r = calc_radius(k, n, b) / 2.0
        theta = 2 * pi * k / phi**2
        rs.append(r)
        ths.append(theta)
    return rs, ths

def calc_radius(k, n, b):
    if k > n - b:
        return 1  # put on the boundary
    else:
        return sqrt(k - 1 / 2) / sqrt(n - (b + 1) / 2)  # apply square root

def fetch_points(radius, method, n = 0):
    if method == "sunflower":
        rs, ts = sunflower(n, 2)
        source = 0
    elif method == "Quad16":
        rs = [0.444,0.230,0.444,0.230,0.444,0.230,0.444,0.230,0.444,
                0.230,0.444,0.230,0.444,0.230,0.444,0.230]
        ts = [0,0.393,0.785,1.178,1.571,1.963,2.356,2.749,3.142,3.534,
                3.927,4.320,4.712,5.105,5.498,5.890]
        source = 0
    elif method == "Cross16":
        xs = [0.047,0.142,0.237,0.332,0,0,0,0,-0.048,-0.144,-0.239,
                -0.334,0,0,0,0]
        ys = [0,0,0,0,0.048,0.143,0.238,0.334,0,0,0,0,-0.048,
                -0.143,-0.239,-0.334]
        source = 1

    if source == 0:
        return [2*radius*r*(cos(t) + 1j * sin(t)) for r, t in zip(rs,ts)]
    elif source == 1:
        return [2*radius*(x + 1j * y) for x, y in zip(xs,ys)]

def veer_effect():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif"
    })
    plt.rc('text.latex', preamble='\\usepackage{braket}')

    fig, axs = plt.subplots(ncols = 4, nrows = 1, figsize=(7,3))
    daoas = [i*pi/180 for i in [5,15,45]]
    alpha_names = ["5", "15", "45"]
    xds = [4,6,8,10]
    yaw = 0*pi/180 
    betas = [0.001,0.5,1,1.5,2,2.5,3,3.5,4]
    betas_refined = linspace(0.0001, 4, 50, endpoint=True)
    markers = ["o","^","s","v","D"]
    delta_handles = []
    delta_labels = []
    circ_handles = []
    anl_handles = []
    num_handles = []
    colors = ["k", "tab:red", "tab:blue", "tab:orange"]
    letters = ["a","b","c","d","e","f","g","h","i","j",
        "k","l","m","n","o","p","q","r","s","t","u",
        "v","w","x","y","z"]
    radius = 124.6 / 2.0 # NREL-5MW
    pss = fetch_points(radius, "sunflower", n = 2000)
    column_index = -1
    ct = 0.8
    ti = 0.05
    delta = 0

    for xd in xds:
        column_index += 1
        ax = axs[column_index]

        # for delta in deltas:
        for daoa in daoas:
            gamma, xi, omega = calc_constants(ct, ti, yaw, xd, daoa)
            id = daoas.index(daoa)
            label = "$\\Delta \\alpha_o = " + alpha_names[id] + "^\circ$"
            ax.text(0.54,0.8-0.07*id,"$\\kappa^{(1)}=" + str(round(calc_kappa(gamma, xi, omega, delta),2)) + "$",
                transform=ax.transAxes,fontsize = 8,
                color = colors[id])
            
            if id < 2:
                curve_circ,  = ax.plot(betas_refined, [solve_new(gamma, beta, xi, omega, delta) for beta in betas_refined],
                    linewidth = 1, color = colors[id],
                    label= label, zorder = 1, linestyle = "--")
            
            curve,  = ax.plot(betas_refined, [solve_new_square(gamma, beta, xi, omega, delta) for beta in betas_refined],
                linewidth = 1, color = colors[id],
                label= label, zorder = 1)
            
            points = ax.scatter(betas, [solve_numerical(pss, radius, gamma, beta, delta, xi, omega) for beta in betas],
                marker=markers[id],
                s = 16,facecolors='white', edgecolor = colors[id],
                zorder = 20, label= label, linewidth = 0.75)

            if column_index == 0:
                delta_handles.append((curve,points))
                delta_labels.append(label)

                anl_handles.append(curve)
                num_handles.append(points)
                if id < 2: circ_handles.append(curve_circ)

        ax.text(0.37,0.9,"(" + letters[column_index] + ") $x/D_o=" + str(xd) + "$",transform=ax.transAxes,fontsize = 10)
        ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
        if column_index == 0:
            ax.set_ylabel("$\\overline{W}^{(1)}/C$", fontsize = 10)
            ax.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
        else:
            ax.yaxis.set_ticklabels([])
        
        ax.set_xticks([0,1,2,3,4])
        ax.set_xlabel("$\\rho/\\sigma$", fontsize = 10)
        ax.tick_params(axis='both', which='major', labelsize=8)

        ax.set_xlim([-0.2,4.2])
        ax.set_ylim([-0.03,0.82])   

    fig.legend(handles = delta_handles, labels = delta_labels, 
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.02, 0.92], ncol=3,
        handlelength = 2.5,
        prop={'size': 10}, frameon=False, columnspacing = 1.0)
    fig.legend(handles = [(circ_handles[0],circ_handles[1])], labels = ["circular"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.57, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    fig.legend(handles = [(anl_handles[0],anl_handles[1],anl_handles[2])], labels = ["rectangular"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.69, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    fig.legend(handles = [(num_handles[0],num_handles[1],num_handles[2])], labels = ["numerical"], 
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.84, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    
    fig.subplots_adjust(bottom=0.135, top=0.92, left=0.075, right=0.975, wspace=0.1, hspace=0.1)
    fig.savefig("veer_effect_rect.png", dpi=500)

def order_effect():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif"
    })
    plt.rc('text.latex', preamble='\\usepackage{braket}')

    fig, axs = plt.subplots(ncols = 4, nrows = 2, figsize=(7,4))
    xds = [4,6,8,10]
    betas = [0.001,0.5,1,1.5,2,2.5,3,3.5,4]
    betas_refined = linspace(0.0001, 4, 50, endpoint=True)
    markers = ["o","^","s","v","D"]
    delta_handles = []
    delta_labels = []
    circ_handles = []
    anl_handles = []
    num_handles = []
    colors = ["k", "tab:red", "tab:blue", "tab:orange"]
    letters = ["a","b","c","d","e","f","g","h","i","j",
        "k","l","m","n","o","p","q","r","s","t","u",
        "v","w","x","y","z"]
    radius = 124.6 / 2.0 # NREL-5MW
    ns = [1,2,3]
    pss = fetch_points(radius, "sunflower", n = 2000)
    column_index = -1
    ct = 0.8
    ti = 0.05
    delta = pi/4 * 0
    daoa = 7*pi/180.0
    yaw = 20*pi/180 

    for xd in xds:
        column_index += 1
        for isol in range(2):
            ax = axs[isol, column_index]

            # for delta in deltas:
            for n in ns:
                gamma, xi, omega = calc_constants(ct, ti, yaw, xd, daoa)
                id = ns.index(n)
                label = "$n = " + str(n) + "$"
                if isol == 0:
                    ax.text(0.06,0.22-0.09*id,"$\\kappa^{(" + str(n) + ")}=" + str(round(calc_kappa(gamma, xi, omega, delta, n=n),2)) + "$",
                        transform=ax.transAxes,fontsize = 8,
                        color = colors[id])

                    curve_circ,  = ax.plot(betas_refined, [solve_new(gamma, beta, xi, omega, delta, n=n) for beta in betas_refined],
                        linewidth = 1, color = colors[id],
                        label= label, zorder = 1, linestyle = "--")
                else:
                    curve,  = ax.plot(betas_refined, [solve_new_square(gamma, beta, xi, omega, delta, n=n) for beta in betas_refined],
                        linewidth = 1, color = colors[id],
                        label= label, zorder = 1)
                    
                points = ax.scatter(betas, [solve_numerical(pss, radius, gamma, beta, delta, xi, omega, n=n) for beta in betas],
                    marker=markers[id],
                    s = 16,facecolors='white', edgecolor = colors[id],
                    zorder = 20, label= label, linewidth = 0.75)

                if column_index == 0:
                    if isol == 0:
                        circ_handles.append(curve_circ)
                    else:
                        delta_handles.append((curve,points))
                        anl_handles.append(curve)
                    delta_labels.append(label)                    
                    num_handles.append(points)
                    

            ax.text(0.37,0.9,"(" + letters[column_index + isol*4] + ") $x/D_o=" + str(xd) + "$",transform=ax.transAxes,fontsize = 10)
            ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
            if column_index == 0:
                if isol == 0:
                    ax.set_ylabel("$\\overline{W}_c^{(n)}/C$", fontsize = 10)
                else: 
                    ax.set_ylabel("$\\overline{W}_r^{(n)}/C$", fontsize = 10)
                ax.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
            else:
                ax.yaxis.set_ticklabels([])
            
            ax.set_xticks([0,1,2,3,4])
            if isol == 1:
                ax.set_xlabel("$\\rho/\\sigma$", fontsize = 10)
            else:
                ax.set_xticklabels([])
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax.set_xlim([-0.2,4.2])
            ax.set_ylim([-0.03,0.82])   

    fig.legend(handles = delta_handles, labels = delta_labels, 
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.04, 0.92], ncol=3,
        handlelength = 2.5,
        prop={'size': 10}, frameon=False, columnspacing = 1.0)
    fig.legend(handles = [(circ_handles[0],circ_handles[1],circ_handles[2])], labels = ["circular"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.55, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    fig.legend(handles = [(anl_handles[0],anl_handles[1],anl_handles[2])], labels = ["rectangular"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.67, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    fig.legend(handles = [(num_handles[0],num_handles[1],num_handles[2])], labels = ["numerical"], 
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.82, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    
    fig.subplots_adjust(bottom=0.1, top=0.92, left=0.075, right=0.975, wspace=0.1, hspace=0.05)
    fig.savefig("order_effect_rect.png", dpi=500)

def splitline(l, char):
    return [x.strip() for x in l.split(char) if x != ""]

def read_windfarm(file):
    lines = (open(file, "r")).readlines()
    yy = [float(splitline(line,"\t")[0]) for line in lines if len(line) > 0]
    xx = [float(splitline(line,"\t")[1]) for line in lines if len(line) > 0]
    x = zeros_like(xx)
    y = zeros_like(xx)
    n = 8
    for i in range(len(xx)):
        columnindex = n - 1 - i%n
        rowindex = int(i/n)
        newindex = columnindex + rowindex*n
        x[newindex] = xx[i]
        y[newindex] = yy[i]
    return x, y

def read_yaw_angles(file):
    yawss = array([[float(x) for x in splitline(line, ",") if len(x) > 0] for line in (open(file,"r")).readlines()])
    n = yawss.shape[0]
    m = yawss.shape[1]
    yaws = []
    for j in range(m):
        for i in range(n):
            ii = n-1-i
            yaws.append(yawss[ii,j])
    return yaws

def read_table(file):
    lines = (open(file)).readlines()[2:]
    return [float(splitline(line, " ")[0]) for line in lines if len(line) > 0], \
           [float(splitline(line, " ")[1]) for line in lines if len(line) > 0], \
           [float(splitline(line, " ")[2]) for line in lines if len(line) > 0]
           
def interpolate_table(u, us, cts):
    if u < us[0]:
        return cts[0]
    elif u > us[-1]:
        return cts[-1]
    else:
        f = interp1d(us, cts)
        return f(u).max()
    
def create_turbine_points(xlocal,ylocal,center):
    points = zeros([len(xlocal),2])
    points[:,0] = center + xlocal
    points[:,1] = ylocal
    return points

def angle_to_vec(angle):
    return array([-sin(angle), -cos(angle)])

def normal_vec(vec):
    return cross(vec, array([0,0,1]))
    
def solve_HornsRev(directory, cordsname, tablename, yawname, sup, is_yawed, anl):
    
    # turbine diameter
    dia = 80
    radius = dia / 2.0

    # wind speed
    u = 8
    tii = 0.077
    
    # wind direction
    wd = 270*pi/180.0
    wdvec = angle_to_vec(wd)
       
    # TI
    # Since this is just to test speed, no need to have a model
    # for turbine-induced turbulence
    # Assume TI is constant throughout the farm
    #ti = 0.077
    #ks = 0.003678 + 0.3837 * ti
    
    # create a square windfarm with num_turbines x num_turbines wind turbines
    xturbines, yturbines = read_windfarm(join(directory, cordsname))
    num_turbines = len(xturbines)
    yaws = read_yaw_angles(join(directory, yawname))
    ut, ctt, _ = read_table(join(directory, tablename))
    
    # create local points on a generic turbine rotor
    xlocal, ylocal = create_rotor_points(dia)

    # initiate wind speed
    speed = [u  for y in yaws]
    speed_num = [u for y in yaws]
    tis = zeros(len(xturbines)) + tii
    
    for it in range(num_turbines):
        xt = xturbines[it]
        yt = yturbines[it]
        points = create_turbine_points(xlocal,ylocal,yt)
        tisq = tis[it]**2

        # loop over upstream turbines
        num_deficit = 1 if sup == "lm" else 0
        anl_deficit = 1 if sup == "lm" else 0
        num_ups = 0
        
        for iup in range(num_turbines):
            if iup == it: continue
            if int(it/8) == int(iup/8): continue # same row (applies for rectangular farms only for normal wind direction)
            xup = xturbines[iup]
            yup = yturbines[iup]
            turbinesvec = array([xt-xup, yt-yup])
            dx = dot(turbinesvec, wdvec)
            if dx <= 0: continue  

            tiup = tis[iup]
            yawup = yaws[iup] * pi/180.0 * is_yawed
            uup = speed[iup] #* cos(yawup)
            ct = interpolate_table(uup, ut, ctt)
            #ti = tis[iup]
            x0 = dia / sqrt(2) * cos(yawup) * (1 + sqrt(1-ct)) / (2.32*tiup+0.154*(1-sqrt(1-ct)))
            th0 = 0.3*yawup/cos(yawup) * (1-sqrt(1-ct*cos(yawup)))
            
            ks = 0.003678 + 0.3837 * tiup
            ky = ks
            kz = ks
            sy = ky*(dx-x0) + dia*cos(yawup)/sqrt(8)
            sz = kz*(dx-x0) + dia/sqrt(8)
            c = 1 - sqrt(1-clip(ct*cos(yawup)/(8*sy*sz/dia**2),0,1))
            if dx < x0:
                deflection = th0 * dx / dia
            else:
                deflection = th0 * x0 / dia + th0/14.7 * sqrt(cos(yawup)/(ky*kz*ct)) * (2.9+1.3*sqrt(1-ct)-ct) * \
                    log( ((1.6+sqrt(ct))*(-sqrt(ct)+1.6*sqrt(8*sy*sz/(dia**2 * cos(yawup)))) ) / ((1.6-sqrt(ct))*(+sqrt(ct)+1.6*sqrt(8*sy*sz/(dia**2 * cos(yawup)))) ))
            deflection *= dia
            #rho0 = norm(turbinesvec - wdvec*dx) + abs(deflection)
            wake_center = yup - deflection
            delta = pi if wake_center < yt else 0
            rho0 = abs(wake_center-yt)
            rho = rho0 / sz
            xi = sqrt(1-(sy/sz)**2)
            
            if rho > 4: continue
            num_ups += 1
            
            # numerical
            dy2 = power(points[:,0]-wake_center, 2)
            dz2 = power(points[:,1], 2)
            deficit =  c * exp(- divide(dy2, 2.0 * sy**2)) * exp(- divide(dz2, 2.0 * sz**2))
            
            if sup == "linear":
                num_deficit += deficit * uup
            elif sup == "rms":
                num_deficit += power((deficit*uup), 2)
            elif sup == "lm":
                num_deficit *= (1-deficit)
            
            # Analytical
            if anl == "circ":
                deficit = solve_new(radius/sz, rho, xi, 0, delta, n=1) * c
            elif anl == "rect":
                deficit = solve_new_square(radius/sz, rho, xi, 0, delta, n=1) * c
            
            if sup == "linear":
                anl_deficit += deficit * uup
            elif sup == "rms":
                anl_deficit += (deficit*uup)**2
            elif sup == "lm":
                anl_deficit *= (1-deficit)
                
            # Add turbulence
            tisq += (crespo_model(uup, axial_induction(ct), dx/dia, rho0/dia))**2

        if num_ups > 0:
            if sup == "lm":
                speed[it] = anl_deficit * u
                speed_num[it] = mean(num_deficit) * u
            elif sup == "linear":
                speed[it] = u - anl_deficit
                speed_num[it] = u - mean(num_deficit)
            elif sup == "rms":
                speed[it] = u - sqrt(anl_deficit)
                speed_num[it] = u - mean(sqrt(num_deficit))
                
            tis[it] = sqrt(tisq)
    
    # Export wind speeds
    f = open(join(directory, "wind_speeds_" + anl + "_" + sup + ("yawed" if is_yawed==1 else "non_yawed") + ".csv"), "w")
    f.write("uanl, unum, \n")
    for u, uu in zip(speed, speed_num): f.write(str(u) + ", " + str(uu) + ",\n")
    f.close()
    
def axial_induction(ct):
    return 0.5*(1-sqrt(clip(1-ct, 0, 1)))

def crespo_model(ws, axial, xd, yd):
    a = 0.8
    b = 0.1
    c = -0.32
    return (0.5 * (axial**a) * (ws**b) * (xd**c)) * int(xd <= 15) * int(yd <= 2)
    
def rotate_turbines(xturbines, yturbines, yaws, width, wd):
    p1s = []
    p2s = []
    for x, y, yaw in zip(xturbines, yturbines, yaws):
        tangent = normal_vec(angle_to_vec((wd - yaw)*pi/180.0))[0:2]
        p1 = array([x,y]) + width/2.0 * tangent
        p2 = array([x,y]) - width/2.0 * tangent
        p1s.append([p1[0], p2[0]])
        p2s.append([p1[1], p2[1]])
    return p1s, p2s   

def row_average(vec):
    return [mean(vec[8*j:8*(j+1)]) for j in range(int(len(vec)/8))] 

def create_figure_no_yaw_rect(directory, yawfile, cordsname, tablename):
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use("Agg")
    
    plt.rcParams.update({"text.usetex": True,"font.family": "sans-serif"})
    
    fig = plt.figure(figsize=(6,3))
    axs = fig.add_gridspec(2,2)
    rows0 = [1,2,3,4,5,6,7,8,9,10]
    xturbines, yturbines = read_windfarm(join(directory, cordsname))
    yaws = read_yaw_angles(join(directory, yawfile))
    p1s, p2s = rotate_turbines(xturbines, yturbines, yaws, 500, 270)
    row_averaged_yaw = row_average(yaws)
    
    # farm layout and yaw angles
    ax = fig.add_subplot(axs[:,0])
    ax.set_xlim([0.5,10.5])
    ax.set_ylim([-1,30])
    ax.plot(rows0, [-y for y in row_averaged_yaw], color = "k")
    ax.scatter(rows0, [-y for y in row_averaged_yaw], marker="o",s = 15, 
                linewidth = 1, facecolors="white", edgecolor = "k", zorder = 10)
    ax.set_xlabel("Row Index", fontsize = 8)
    ax.set_ylabel("$\\gamma~(^\circ)$", fontsize = 8)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_xticks(rows0)
    ax.set_xticklabels(rows0)
    ax.text(0.08, 0.93, "(a)", transform=ax.transAxes, fontsize = 9)
    
    xc = mean(xturbines)
    yc = mean(yturbines)
    rs = max([sqrt((xx-xc)**2 + (yy-yc)**2) for xx, yy in zip(xturbines, yturbines)])*1.15

    ths = linspace(0,2*pi,50,endpoint=True)
    circle_x = [rs*cos(th)+xc for th in ths]
    circle_y = [rs*sin(th)+yc for th in ths]
    turbine_colors = ["k" for x in xturbines]
    t1 = turbine_colors.copy()
    for kc in range(8):
        t1[kc] = "tab:green"

    bounds = [-1.91,12.3,17,17]
    axi = ax.inset_axes(bounds=bounds, zorder = 100, transform=ax.transData)
    
    for p1, p2, c in zip(p1s, p2s, t1):
        axi.plot(p1, p2, color= c, linewidth = 0.5, zorder = 5)
    
    axi.scatter(xturbines,yturbines, marker="o",s = 4, 
                linewidth = 0.5, facecolors="white", edgecolor = t1, zorder = 10)
    axi.set_aspect("equal")
    axi.plot(circle_x, circle_y, color = "k", linestyle="--", linewidth = 0.75)
    axi.get_xaxis().set_visible(False)
    axi.get_yaxis().set_visible(False)
    axi.spines['left'].set_visible(False)
    axi.spines['bottom'].set_visible(False)
    axi.spines['right'].set_visible(False)
    axi.spines['top'].set_visible(False)
    axi.arrow(xc-1.1*rs,yc, 750, 0, head_width = 200, facecolor="k")
    
    ## power generation
    names = ["Linear", "RMS", "Product"]
    colors = ["k", "tab:red", "tab:blue"]
    sup_handles = []
    labels = []
    circ_handles = []
    rect_handles = []
    num_handles = []
    markers = ["o","^","s"]
    for isol, style in zip(range(2),["--","-"]):
        if isol == 0:
            files = ["wind_speeds_circ_linearyawed.csv",
                     "wind_speeds_circ_rmsyawed.csv",
                     "wind_speeds_circ_lmyawed.csv"]
        else:
            files = ["wind_speeds_rect_linearyawed.csv",
                     "wind_speeds_rect_rmsyawed.csv",
                     "wind_speeds_rect_lmyawed.csv"]
        ax = fig.add_subplot(axs[isol,1])
        ax.set_ylabel("$P/P_1$", fontsize = 8)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.set_xticks(rows0)
        ax.set_xlim([0.5,10.5])
        ax.set_ylim([0.6,1.02])
        ut, _, pt = read_table(join(directory, tablename))
        yaws_avg = row_average(yaws)
        for file, name, c, mark in zip(files, names, colors, markers):
            uanl, unum = read_speed_file(join(directory, file))
            uanl_avg = row_average([ua for ua in uanl])
            unum_avg = row_average([ua for ua in unum])
            
            panl = [interpolate_table(u, ut, pt)* (cos(1*y*pi/180.0))**1.8 for u, y in zip(uanl_avg,yaws_avg)]
            pnum = [interpolate_table(u, ut, pt)* (cos(1*y*pi/180.0))**1.8 for u, y in zip(unum_avg,yaws_avg)]
            curve,  = ax.plot(rows0, [p/panl[0] for p in panl],
                            color = c, linewidth = 1, linestyle = style, 
                            label=name, zorder = 10 - names.index(name))
        
            points = ax.scatter(rows0, [p/pnum[0] for p in pnum], marker=mark,s = 15, 
                    linewidth = 1, facecolors="white", edgecolor = c, 
                    zorder = 50 - names.index(name))

            if isol == 0: 
                title = "(b) circular"
                ax.set_xticklabels([])
                circ_handles.append(curve)
            else:
                sup_handles.append((curve, points))
                labels.append(name)
                title = "(c) rectangular"
                ax.set_xlabel("Row Index", fontsize = 8)
                ax.set_xticklabels(rows0)
                rect_handles.append(curve)
                num_handles.append(points)
                
        ax.text(0.65, 0.86, title, transform=ax.transAxes, fontsize = 9)

    fig.legend(handles = sup_handles, labels = labels,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.04, 0.92], ncol=3,
        handlelength = 2.5,
        prop={'size': 9}, frameon=False, columnspacing = 1.0)
    fig.legend(handles = [(circ_handles[0],circ_handles[1],circ_handles[2])], labels = ["circular"],
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.5,
        loc = [0.53, 0.92], ncol=4, prop={'size': 9}, frameon=False)
    fig.legend(handles = [(rect_handles[0],rect_handles[1],rect_handles[2])], labels = ["rectangular"],
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.66, 0.92], ncol=4, prop={'size': 9}, frameon=False)
    fig.legend(handles = [(num_handles[0],num_handles[1],num_handles[2])], labels = ["numerical"],
        handler_map={tuple: HandlerTuple(ndivide=None)},
        handlelength = 2,
        loc = [0.82, 0.92], ncol=4, prop={'size': 9}, frameon=False)
    fig.subplots_adjust(bottom=0.12, top=0.92, left=0.075, right=0.975, wspace=0.22, hspace=0.1)
    fig.savefig(join(directory, "Horns_yawed_270_rect.png"), dpi=500)    

def uncertainity():
    cases = [["Rect", 0],
            #["Circ", 0],
            ["Quad16", 16],
            ["Cross16", 16],
            ["sunflower", 16],
            ["sunflower", 100],
            ["sunflower", 500],
            ["sunflower", 1000]]
    
    radius = 124.6 / 2.0 # NREL-5MW
    betas = [0.001,0.5,1,1.5,2,2.5,3,3.5,4]
    #daoas = [i*pi/180 for i in [0.001,1,2,3,4,5,6,7]]
    daoas = [i*pi/180 for i in [10,15,20,25,30,35,40,45]]
    yaws = [i*pi/180 for i in [0.01,10,20,30]]
    ns = [1,2,3]
    deltas = [0,pi/4,pi/2,3*pi/4]
    xds = [4,6,8,10]
    
    pss_ref = fetch_points(radius, "sunflower", n = 2000)   
    ct = 0.8
    ti = 0.05
    rmse = [0 for c in cases]
    count = 0

    for yaw in yaws:
        for daoa in daoas:
            for n in ns:
                for xd in xds:
                    for delta in deltas:
                        gamma, xi, omega = calc_constants(ct, ti, yaw, xd, daoa)
                        ref_sol = [solve_numerical(pss_ref, radius, gamma, beta, delta, xi, omega, n=n) for beta in betas]

                        icase = -1
                        for case in cases:
                            icase += 1
                            if case[0] == "Rect":
                                sol = [solve_new_square(gamma, beta, xi, omega, delta, n=n) for beta in betas]
                            elif case[0] == "Circ":
                                sol = [solve_new(gamma, beta, xi, omega, delta, n=n) for beta in betas]
                            else:
                                pss = fetch_points(radius, case[0], n = case[1])   
                                sol = [solve_numerical(pss, radius, gamma, beta, delta, xi, omega, n=n) for beta in betas]

                            rmse[icase] += calc_rmse(ref_sol, sol)
                            count += len(betas)
    icase = -1
    for case in cases:
        icase += 1
        print(case[0] + " " + str(case[1]) + ": " + str(1000*sqrt(rmse[icase] / count)))

def calc_rmse(ref_sol, sol):
    q = array(ref_sol) - array(sol)
    return dot(q,q)

def plot_points_2():
    cases = [["Quad16", 16, "Q16"],
            ["Cross16", 16, "C16"],
            ["sunflower", 16, "S16"],
            ["sunflower", 100, "S100"],
            ["sunflower", 500, "S500"],
            ["sunflower", 1000, "S1000"]]
    
    letters = ["a","b","c","d","e","f","g","h","i","j",
        "k","l","m","n","o","p","q","r","s","t","u",
        "v","w","x","y","z"]

    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif"
    })
    plt.rc('text.latex', preamble='\\usepackage{braket}')

    radius = 0.5
    fig, axs = plt.subplots(ncols = 3, nrows = 2, figsize=(6,4))
    axs = axs.flatten()
    thetas = linspace(0,2*pi,100,endpoint=True)
    sizes = [3,3,3,3,1,1]

    for icase in range(len(cases)):
        case = cases[icase]
        ax = axs[icase]
        pss = fetch_points(radius, case[0], n = case[1])
        ax.scatter([p.real for p in pss], [p.imag for p in pss], marker="o",
                s = sizes[icase], edgecolor = "tab:red", zorder = 20,
                label= case[2], linewidth = 0.75, facecolor = "tab:red")
        ax.plot([radius*cos(t) for t in thetas], [radius*sin(t) for t in thetas],
        linestyle = "--", color = "k", linewidth = 0.75)
        ax.text(0.8, 0.9, "(" + letters[icase] + ") " +  case[2], transform=ax.transAxes,fontsize = 8)

        ax.set_xlim([-0.55,0.55])
        ax.set_xticks([-0.5,0,0.5])
        ax.set_yticks([-0.5,0,0.5])
        ax.set_ylim([-0.55,0.55])
        if icase in [0,3]:
            ax.set_ylabel("$z/D$", fontsize = 9)
        else:
            ax.set_yticklabels([])
        if icase > 2:
            ax.set_xlabel("$y/D$", fontsize = 9)
        else:
            ax.set_xticklabels([])
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.set_aspect("equal")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    
    fig.subplots_adjust(bottom=0.13, top=0.95, left=0.075, right=0.975, wspace=0.05, hspace=0.1)
    fig.savefig("points.png", dpi=500)

def validate_hub_height():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif"
    })
    plt.rc('text.latex', preamble='\\usepackage{braket}')
    
    fig, axs = plt.subplots(ncols = 4, nrows = 1, figsize=(7,3))
    daoa = 7*pi/180
    xds = [4,6,8,10]
    yaw = 20*pi/180 
    betas = [0.001,0.5,1,1.5,2,2.5,3,3.5,4]
    betas_refined = linspace(0.0001, 4, 50, endpoint=True)
    deltas = [0, pi/4, 3*pi/4]
    markers = ["o","^","s","v","D"]
    delta_handles = []
    delta_labels = []
    anl_handles = []
    num_handles = []
    colors = ["k", "tab:red", "tab:blue", "tab:orange"]
    letters = ["a","b","c","d","e","f","g","h","i","j",
        "k","l","m","n","o","p","q","r","s","t","u",
        "v","w","x","y","z"]
    delta_names = ["0", "\\pi/4", "3\\pi/4"]
    radius = 124.6 / 2.0 # NREL-5MW
    column_index = -1
    ct = 0.8
    ti = 0.05

    for xd in xds:
        gamma, xi, omega = calc_constants(ct, ti, yaw, xd, daoa)
        column_index += 1
        ax = axs[column_index]

        for delta in deltas:
            idelta = deltas.index(delta)
            label = "$\\delta=" + delta_names[idelta] + "$"
            
            curve,  = ax.plot(betas_refined, [solve_new(gamma, beta, xi, omega, delta) for beta in betas_refined],
                linewidth = 1, color = colors[idelta],
                label= label, zorder = 1)
            
            points, = ax.plot(betas_refined, [solve_numerical([0], radius, gamma, beta, delta, xi, omega) for beta in betas_refined],
                linewidth = 0.75, color = colors[idelta],
                label= label, zorder = 1, linestyle="--")
            
            if column_index == 0:
                delta_handles.append((curve,points))
                delta_labels.append(label)

                anl_handles.append(curve)
                num_handles.append(points)
            ax.text(0.37,0.9,"(" + letters[column_index] + ") $x/D_o=" + str(xd) + "$",transform=ax.transAxes,fontsize = 10)
            ax.text(0.54,0.8,"$\\xi=" + str(round(xi,2)) + "$",transform=ax.transAxes,fontsize = 8)
            ax.text(0.54,0.73,"$\\omega=" + str(round(omega,2)) + "$",transform=ax.transAxes,fontsize = 8)
            ax.text(0.54,0.66,"$R/\\sigma=" + str(round(gamma,1)) + "$",transform=ax.transAxes,fontsize = 8)
            ax.text(0.54,0.59,"$\\kappa^{(1)}=" + str(round(calc_kappa(gamma, xi, omega, delta),2)) + "$",transform=ax.transAxes,fontsize = 8)
            ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
            if column_index == 0:
                ax.set_ylabel("$\\braket{\\overline{W}_c^{(1)},\\hat{W}}/C$", fontsize = 9)
                ax.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
            else:
                ax.yaxis.set_ticklabels([])
            
            ax.set_xticks([0,1,2,3,4])
            ax.set_xlabel("$\\rho/\\sigma$", fontsize = 10)
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax.set_xlim([-0.2,4.2])
            ax.set_ylim([-0.03,1.02])   

    fig.legend(handles = delta_handles, labels = delta_labels, 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        loc = [0.06, 0.92], ncol=3,
        handlelength = 2,
        prop={'size': 10}, frameon=False)
    fig.legend(handles = [(anl_handles[0],anl_handles[1],anl_handles[2])], labels = ["rotor-averaged"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.66, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    fig.legend(handles = [(num_handles[0],num_handles[1],num_handles[2])], labels = ["nacelle"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        loc = [0.85, 0.92], 
        handlelength = 1.4,
        ncol=4, prop={'size': 10}, frameon=False)
    
    fig.subplots_adjust(bottom=0.135, top=0.92, left=0.075, right=0.975, wspace=0.1, hspace=0.1)
    fig.savefig("validate_circle_single_nacelle.png", dpi=500)

def yaw_effect():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif"
    })
    plt.rc('text.latex', preamble='\\usepackage{braket}')

    fig, axs = plt.subplots(ncols = 4, nrows = 2, figsize=(7,4))
    xds = [4,6,8,10]
    yaws = [q * pi/180 for q in [0.0001,15,30]]
    yaw_names = ["0","15","30"]
    betas = [0.001,0.5,1,1.5,2,2.5,3,3.5,4]
    betas_refined = linspace(0.0001, 4, 50, endpoint=True)
    markers = ["o","^","s","v","D"]
    delta_handles = []
    delta_labels = []
    circ_handles = []
    anl_handles = []
    num_handles = []
    colors = ["k", "tab:red", "tab:blue", "tab:orange"]
    letters = ["a","b","c","d","e","f","g","h","i","j",
        "k","l","m","n","o","p","q","r","s","t","u",
        "v","w","x","y","z"]
    radius = 124.6 / 2.0 # NREL-5MW
    pss = fetch_points(radius, "sunflower", n = 2000)
    column_index = -1

    ct = 0.8
    ti = 0.05
    delta = 0
    daoa = 0*pi/180.0

    for xd in xds:
        column_index += 1

        for isol in range(2):
            ax = axs[isol, column_index]

            # for delta in deltas:
            for yaw in yaws:
                gamma, xi, omega = calc_constants(ct, ti, yaw, xd, daoa)
                id = yaws.index(yaw)
                label = "$\\gamma_o = " + yaw_names[id] + "^\circ$"
                ax.text(0.63,0.75-0.09*id,"$\\xi=" + str(round(xi,2)) + "$",
                transform=ax.transAxes,fontsize = 8,
                color = colors[id])
            
                if isol == 0:
                    curve_circ,  = ax.plot(betas_refined, [solve_new(gamma, beta, xi, omega, delta) for beta in betas_refined],
                        linewidth = 1, color = colors[id],
                        label= label, zorder = 1, linestyle = "--")
                else:
                    curve,  = ax.plot(betas_refined, [solve_new_square(gamma, beta, xi, omega, delta) for beta in betas_refined],
                        linewidth = 1, color = colors[id],
                        label= label, zorder = 1)
                
                points = ax.scatter(betas, [solve_numerical(pss, radius, gamma, beta, delta, xi, omega) for beta in betas],
                    marker=markers[id],
                    s = 16,facecolors='white', edgecolor = colors[id],
                    zorder = 20, label= label, linewidth = 0.75)

                if column_index == 0:
                    if isol == 0:
                        circ_handles.append(curve_circ)
                    else:
                        anl_handles.append(curve)
                        delta_handles.append((curve,points))
                    delta_labels.append(label)
                    num_handles.append(points)

            ax.text(0.37,0.9,"(" + letters[column_index + isol*4] + ") $x/D_o=" + str(xd) + "$",transform=ax.transAxes,fontsize = 10)
            ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
            if column_index == 0:
                if isol == 0:
                    ax.set_ylabel("$\\overline{W}_c^{(1)}/C$", fontsize = 10)
                else:
                    ax.set_ylabel("$\\overline{W}_r^{(1)}/C$", fontsize = 10)
                ax.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
            else:
                ax.yaxis.set_ticklabels([])
            
            ax.set_xticks([0,1,2,3,4])
            if isol == 1:
                ax.set_xlabel("$\\rho/\\sigma$", fontsize = 10)
            else:
                ax.set_xticklabels([])
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax.set_xlim([-0.2,4.2])
            ax.set_ylim([-0.03,0.82])   

    fig.legend(handles = delta_handles, labels = delta_labels, 
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.02, 0.92], ncol=3,
        handlelength = 2.5,
        prop={'size': 10}, frameon=False, columnspacing = 1.0)
    fig.legend(handles = [(circ_handles[0],circ_handles[1],circ_handles[2])], labels = ["circular"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.57, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    fig.legend(handles = [(anl_handles[0],anl_handles[1],anl_handles[2])], labels = ["rectangular"], 
        handler_map={tuple: HandlerTupleVertical(ndivide=None)},
        handlelength = 1.2,
        loc = [0.69, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    fig.legend(handles = [(num_handles[0],num_handles[1],num_handles[2])], labels = ["numerical"], 
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc = [0.84, 0.92], ncol=4, prop={'size': 10}, frameon=False)
    
    fig.subplots_adjust(bottom=0.12, top=0.92, left=0.075, right=0.975, wspace=0.1, hspace=0.1)
    fig.savefig("yaw_effect.png", dpi=500)
    
def create_figure_no_yaw(directory, linearfile, rmsfile, lmfile, yawfile,
                  cordsname, tablename):
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use("Agg")
    
    plt.rcParams.update({"text.usetex": True,"font.family": "sans-serif"})
    
    fig, axs = plt.subplots(ncols = 2, nrows = 1, figsize=(6,3))    
    rows0 = [1,2,3,4,5,6,7,8,9,10]
    xturbines, yturbines = read_windfarm(join(directory, cordsname))
    yaws = read_yaw_angles(join(directory, yawfile))
    p1s, p2s = rotate_turbines(xturbines, yturbines, yaws, 500, 270)
    row_averaged_yaw = row_average(yaws)
    
    # farm layout and yaw angles
    ax = axs[0]
    ax.set_xlim([0.5,10.5])
    ax.set_ylim([-1,30])
    ax.plot(rows0, [-y for y in row_averaged_yaw], color = "k")
    ax.scatter(rows0, [-y for y in row_averaged_yaw], marker="o",s = 15, 
                linewidth = 1, facecolors="white", edgecolor = "k", zorder = 10)
    ax.set_xlabel("Row Index", fontsize = 8)
    ax.set_ylabel("$\\gamma~(^\circ)$", fontsize = 8)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_xticks(rows0)
    ax.set_xticklabels(rows0)
    ax.text(0.08, 0.93, "(a)", transform=ax.transAxes, fontsize = 8)
    
    xc = mean(xturbines)
    yc = mean(yturbines)
    rs = max([sqrt((xx-xc)**2 + (yy-yc)**2) for xx, yy in zip(xturbines, yturbines)])*1.15

    ths = linspace(0,2*pi,50,endpoint=True)
    circle_x = [rs*cos(th)+xc for th in ths]
    circle_y = [rs*sin(th)+yc for th in ths]
    turbine_colors = ["k" for x in xturbines]
    t1 = turbine_colors.copy()
    for kc in range(8):
        t1[kc] = "tab:green"

    bounds = [-1.91,12.3,17,17]
    axi = ax.inset_axes(bounds=bounds, zorder = 100, transform=ax.transData)
    
    for p1, p2, c in zip(p1s, p2s, t1):
        axi.plot(p1, p2, color= c, linewidth = 0.5, zorder = 5)
    
    axi.scatter(xturbines,yturbines, marker="o",s = 4, 
                linewidth = 0.5, facecolors="white", edgecolor = t1, zorder = 10)
    axi.set_aspect("equal")
    axi.plot(circle_x, circle_y, color = "k", linestyle="--", linewidth = 0.75)
    axi.get_xaxis().set_visible(False)
    axi.get_yaxis().set_visible(False)
    axi.spines['left'].set_visible(False)
    axi.spines['bottom'].set_visible(False)
    axi.spines['right'].set_visible(False)
    axi.spines['top'].set_visible(False)
    axi.arrow(xc-1.1*rs,yc, 750, 0, head_width = 200, facecolor="k")
    
    ## power generation
    files = [linearfile, rmsfile, lmfile]
    names = ["Linear", "RMS", "Product"]
    colors = ["k", "tab:red", "tab:blue"]
    ax = axs[1]
    ax.set_ylabel("$P/P_1$", fontsize = 8)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_xticks(rows0)
    ax.set_xlabel("Row Index", fontsize = 8)
    ax.set_xticklabels(rows0)
    ax.set_xlim([0.5,10.5])
    ax.set_ylim([0.6,1.02])
    
    ut, _, pt = read_table(join(directory, tablename))
    yaws_avg = row_average(yaws)
    handles = []
    for file, name, c, item in zip(files, names, colors, ["b","c","d"]):
        uanl, unum = read_speed_file(join(directory, file + "yawed" + ".csv"))
        uanl_avg = row_average([ua for ua in uanl])
        unum_avg = row_average([ua for ua in unum])
        
        panl = [interpolate_table(u, ut, pt)* (cos(1*y*pi/180.0))**1.8 for u, y in zip(uanl_avg,yaws_avg)]
        pnum = [interpolate_table(u, ut, pt)* (cos(1*y*pi/180.0))**1.8 for u, y in zip(unum_avg,yaws_avg)]
        curve,  = ax.plot(rows0, [p/panl[0] for p in panl],
                        color = c, linewidth = 1, linestyle = "-", 
                        label=name, zorder = 10 - names.index(name))
        
        handles.append(curve)
        ax.scatter(rows0, [p/pnum[0] for p in pnum], marker="o",s = 15, 
                linewidth = 1, facecolors="white", edgecolor = c, 
                zorder = 50 - names.index(name))
    ax.legend(handles = handles, loc = [0.7, 0.8], frameon = False, ncol=1, prop={'size': 7})
    ax.text(0.08, 0.93, "(b)", transform=ax.transAxes, fontsize = 8)
    fig.subplots_adjust(bottom=0.12, top=0.95, left=0.075, right=0.975, wspace=0.22, hspace=0.2)
    fig.savefig(join(directory, "Horns_yawed_270_simple.pdf"), dpi=500)    
    
def create_figure(directory, linearfile, rmsfile, lmfile, yawfile,
                  cordsname, tablename):
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use("Agg")
    
    plt.rcParams.update({"text.usetex": True,"font.family": "sans-serif"})
    
    #fig, axs = plt.subplots(ncols = 2, nrows = 1, figsize=(6,3))
    fig = plt.figure(figsize=(6,3))
    axs = fig.add_gridspec(3,2)
    rows0 = [1,2,3,4,5,6,7,8, 9, 10]
    xturbines, yturbines = read_windfarm(join(directory, cordsname))
    yaws = read_yaw_angles(join(directory, yawfile))
    p1s, p2s = rotate_turbines(xturbines, yturbines, yaws, 500, 270)
    row_averaged_yaw = row_average(yaws)
    
    # farm layout and yaw angles
    #ax = axs[0]
    ax = fig.add_subplot(axs[:,0])
    ax.set_xlim([0.5,10.5])
    ax.set_ylim([-1,30])
    ax.plot(rows0, [-y for y in row_averaged_yaw], color = "k")
    ax.scatter(rows0, [-y for y in row_averaged_yaw], marker="o",s = 15, 
                linewidth = 1, facecolors="white", edgecolor = "k", zorder = 10)
    ax.set_xlabel("Row Index", fontsize = 8)
    ax.set_ylabel("$\\gamma~(^\circ)$", fontsize = 8)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_xticks(rows0)
    ax.set_xticklabels(rows0)
    ax.text(0.08, 0.95, "(a)", transform=ax.transAxes, fontsize = 8)
    
    xc = mean(xturbines)
    yc = mean(yturbines)
    rs = max([sqrt((xx-xc)**2 + (yy-yc)**2) for xx, yy in zip(xturbines, yturbines)])*1.15

    ths = linspace(0,2*pi,50,endpoint=True)
    circle_x = [rs*cos(th)+xc for th in ths]
    circle_y = [rs*sin(th)+yc for th in ths]
    turbine_colors = ["k" for x in xturbines]
    t1 = turbine_colors.copy()
    for kc in range(8):
        t1[kc] = "tab:green"

    bounds = [-1.91,12.3,17,17]
    axi = ax.inset_axes(bounds=bounds, zorder = 100, transform=ax.transData)
    
    for p1, p2, c in zip(p1s, p2s, t1):
        axi.plot(p1, p2, color= c, linewidth = 0.5, zorder = 5)
    
    axi.scatter(xturbines,yturbines, marker="o",s = 4, 
                linewidth = 0.5, facecolors="white", edgecolor = t1, zorder = 10)
    axi.set_aspect("equal")
    axi.plot(circle_x, circle_y, color = "k", linestyle="--", linewidth = 0.75)
    axi.get_xaxis().set_visible(False)
    axi.get_yaxis().set_visible(False)
    axi.spines['left'].set_visible(False)
    axi.spines['bottom'].set_visible(False)
    axi.spines['right'].set_visible(False)
    axi.spines['top'].set_visible(False)
    axi.arrow(xc-1.1*rs,yc, 750, 0, head_width = 200, facecolor="k")
    
    ## power generation
    files = [linearfile, rmsfile, lmfile]
    names = ["Linear", "RMS", "Product"]
    colors = ["k", "tab:red", "tab:blue"]
    
    axs2 = []
    for i in range(len(names)):
        ax = fig.add_subplot(axs[i,1])
        axs2.append(ax)
        ax.set_ylabel("$P~(\\textrm{MW})$", fontsize = 8)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.set_xticks(rows0)
        if i == 2:
            ax.set_xlabel("Row Index", fontsize = 8)
            ax.set_xticklabels(rows0)
        else:
            ax.set_xticklabels([])
        ax.set_xlim([0.5,10.5])
        ax.set_ylim([0.3,0.8])
    
    namess = {"yawed":"yaw", "non_yawed":"no yaw"}
    ut, _, pt = read_table(join(directory, tablename))
    yaws_avg = row_average(yaws)
    for file, name, c, ax, item in zip(files, names, colors, axs2, ["b","c","d"]):
        handles = []
        for yname, style, g, m in zip(["yawed", "non_yawed"], ["-", "--"], [1,0], ["o", "s"]):
            uanl, unum = read_speed_file(join(directory, file + yname + ".csv"))
            uanl_avg = row_average([ua for ua in uanl])
            unum_avg = row_average([ua for ua in unum])
            
            panl = [interpolate_table(u, ut, pt)* (cos(g*y*pi/180.0))**1.8 for u, y in zip(uanl_avg,yaws_avg)]
            pnum = [interpolate_table(u, ut, pt)* (cos(g*y*pi/180.0))**1.8 for u, y in zip(unum_avg,yaws_avg)]
            curve,  = ax.plot(rows0, [p/1000 for p in panl],
                            color = c, linewidth = 1, linestyle = style, 
                            label=namess[yname], zorder = 10 - names.index(name))
           
            handles.append(curve)
            if style=="-": ax.scatter(rows0, [p/1000 for p in pnum], marker=m,s = 15, 
                    linewidth = 1, facecolors="white", edgecolor = c, 
                    zorder = 50 - names.index(name))
            ax.legend(handles = handles, loc = [0.7, 0.52], frameon = False, ncol=1, prop={'size': 7})
            ax.text(0.09, 0.8, "(" + item + ") " + name, transform=ax.transAxes, fontsize = 7)
    
    fig.subplots_adjust(bottom=0.12, top=0.95, left=0.075, right=0.975, wspace=0.22, hspace=0.2)
    
    fig.savefig(join(directory, "Horns_yawed_270.png"), dpi=500)
    
def use_yaw():
    return 0.000001
    #return 1
    
def read_speed_file(file):
    lines = (open(file, "r")).readlines()[1:]
    return [float(splitline(line,",")[0]) for line in lines if len(line) > 0], \
        [float(splitline(line,",")[1]) for line in lines if len(line) > 0],

def test_speed_nonsym(nx, ny, sx, sy, epoch, nps_in = 50, flag = 0):
    # turbine diameter
    d = 220
    radius = d / 2.0

    # wind speed
    u = 10
    
    # TI
    # Since this is just to test speed, no need to have a model
    # for turbine-induced turbulence
    # Assume TI is constant throughout the farm
    ti = 6/100.0
    k = 0.003678 + 0.3837 * ti
    
    xi = 0.2
    omega = 0.1
    
    # assume all have the same hub-height
    delta = 0

    o2 = omega * omega
    xi2 = xi * xi

    cdelta = cos(delta)
    sdelta = sin(delta)
    tdelta = sdelta / cdelta

    one_m_xi2 = 1 - xi2
    half_o2_m_xi2 = (o2 - xi2) / 2.0

    sp_sq = one_m_xi2 / (1 + half_o2_m_xi2)
    sns_sq = one_m_xi2 /  sqrt(half_o2_m_xi2 * half_o2_m_xi2 + o2)  # **(0.5)
    phi_ns =  arctan2(2 * omega, xi2 - o2)
    phi_s =  arctan2(omega + one_m_xi2 * tdelta / (omega * tdelta + 1), 1)
    ss_sq = one_m_xi2 *  cos(phi_s) / (cdelta + omega * sdelta)
    ss_sq2 = ss_sq * ss_sq

    phi = 2 * phi_s - phi_ns
    phicos =  cos(phi)
    phisin =  sin(phi)

    # create a square windfarm with num_turbines x num_turbines wind turbines
    xturbines, yturbines = create_farm(nx, ny, sx, sy, d)
    
    # create local points on a generic turbine rotor
    nums = []
    # npss = [16, 50, 100, 500, 1000, 1500, 2000]
    npss = [nps_in]
    cases = ["flower" for n in npss]
    for label, nps in zip(cases, npss):
        xlocal, ylocal, weights = create_rotor_points(d, label, nps)
        nums.append([xlocal, ylocal, weights])  
    
    # initiate wind speed
    speed =  zeros(len(xturbines)) + u
    
    index = -1
    time_anl = 0
    times_num = [0 for kk in cases]
    
    for i in range(nx):
        for j in range(ny):
            index += 1
            xt = xturbines[index]
            yt = yturbines[index]

            # loop over upstream turbines
            upindex = -1
            deficit_anl = 0
            defs_num = [0 for kk in range(len(cases))]
            num_ups = 0
            for iup in range(i):
                for jup in range(ny):
                    upindex += 1
                    xup = xturbines[upindex]
                    yup = yturbines[upindex]
                    if xup > xt: continue
                    dx = xt - xup
                    dy = abs(yt - yup)

                    ct = interpolate_ct(speed[upindex])
                    epsilon = 0.2 *  sqrt(0.5*((1+ sqrt(1-ct))/ sqrt(1-ct)))
                    sigma = k*dx + epsilon * d
                    c = 1- sqrt(1- clip(ct/(8*(sigma/d)**2),0,1))

                    if dy / sigma > 5: continue
                    num_ups += 1
                    sigmasqr = sigma**2
                    two_sigmasqr = 2 * sigmasqr
                    
                    # All previous lines are shared by the two approaches
                    # Hence, they are not included in the time measurement
                    
                    ## Analytical approach
                    if flag == 0:
                        ## Circular disk
                        local_deficit = 0
                        start_anl = time.process_time()
                        for t in range(epoch):
                            
                            gamma = radius/sigma
                            g2 = gamma**2
                            
                            rho = dy/sigma # as delta was assumed 0
                            if rho == 0: rho = 0.0001

                            cdelta =  cos(delta)
                            sdelta =  sin(delta)
                            tdelta = sdelta / cdelta
                            one_m_xi2 = 1-xi2
                            half_o2_m_xi2 = (o2 - xi2) / 2.0

                            sp_sq = one_m_xi2 / ( 1 + half_o2_m_xi2 )
                            sns_sq = one_m_xi2 / (half_o2_m_xi2**2 + o2)**(0.5)
                            phi_ns =   atan2(2*omega, xi2 - o2)
                            phi_s =   atan2(omega + one_m_xi2 * tdelta / (omega*tdelta + 1), 1)
                            ss_sq = one_m_xi2 *  cos(phi_s) / (cdelta + omega * sdelta)

                            aux1 = sp_sq * rho
                            aux2 = rho/ss_sq
                            chi2 = aux1*aux2*sp_sq/4.0/sns_sq/ss_sq
                            phi = 2*phi_s - phi_ns
                            pns =  exp(-chi2 *  cos(phi)) *  cos(chi2 *  sin(phi)) - 1
                            lam_rho = gamma * ss_sq / ( aux1 )
                            kappa = gamma*aux2
                            ii0 = i0(kappa)
                            ii1 = i1(kappa)
                            # ii2 = iv(2,kappa)
                            ii2 = ii0 - 2/kappa * ii1
                            R_ss2 = 0.5 * g2/ sp_sq
                            exp_R_ss2 =  exp(-R_ss2)
                            mu0 = sp_sq/g2 * exp_R_ss2 * generic_Phi(R_ss2, ii0, kappa*ii1,(aux2**2 *sp_sq), 20)
                            G = mu0*(1+2*pns) -  exp(-R_ss2) / R_ss2 * pns * (lam_rho * ii1 + lam_rho*lam_rho * ii2)
                            shat_sq = (1/sp_sq +  cos(2*delta-phi_ns)/sns_sq)
                            local_deficit += speed[upindex] * 2 * c *  exp(-rho*rho/2*shat_sq)*G
                        deficit_anl += local_deficit/epoch
                        end_anl = time.process_time()
                        time_anl += end_anl - start_anl
                    
                    elif flag == 1:
                        ## rectangular
                        local_deficit = 0
                        start_anl = time.process_time()
                        for zz in range(epoch):
                            
                            gamma = radius/sigma
                            rho = dy/sigma                            
                            gz = gamma * 0.9
                            gy = gamma * 0.9
                            a = omega/  sqrt(1-xi*xi)
                            rhosd = rho*  sin(delta)
                            rhocd = rho*  cos(delta)
                            qq =   sqrt(1-xi*xi)
                            t = + gen_owen_t(rhosd - gz, a, (rhocd + gy)/ qq) \
                                - gen_owen_t(rhosd + gz, a, (rhocd + gy)/ qq) \
                                - gen_owen_t(rhosd - gz, a, (rhocd - gy)/ qq) \
                                + gen_owen_t(rhosd + gz, a, (rhocd - gy)/ qq)
                            local_deficit += speed[upindex]  * c * (  pi* qq / (2*gy*gz) * t)
                    
                        deficit_anl += local_deficit/epoch
                        end_anl = time.process_time()
                        time_anl += end_anl - start_anl

                    ## Numerical approaches
                    jhh = -1
                    for numcase in nums:
                        jhh += 1
                        deficit = 0
                        start_num = time.process_time()
                        for t in range(epoch):
                            ys = yt -  array(numcase[0]) - yup
                            zs =  array(numcase[1])
                            xi2 = xi ** 2
                            one_m_xi2 = 1 - xi2
                            sysqr = two_sigmasqr*one_m_xi2
                            deficit += speed[upindex] * c *  sum(numcase[2] *  exp(-  divide( power(ys+ multiply(zs,omega), 2), sysqr)) *  exp(-  divide( power(zs, 2), two_sigmasqr)))
                        deficit /= epoch
                        end_num = time.process_time()
                        defs_num[jhh] = defs_num[jhh] + deficit
                        times_num[jhh] = times_num[jhh] + end_num - start_num

            if num_ups > 0: speed[index] = (u - deficit_anl / num_ups)

    #line = float_to_string(time_anl) + " "
    #for d in time_anl: line += float_to_string(d) + " " 
    print('anl',time_anl)
    
    #line = float_to_string(time_anl) + " " 
    #for t in times_num: line += float_to_string(t) + " "
    print('num',times_num)

def test_speed_sym(nx, ny, sx, sy, epoch):
    # turbine diameter
    d = 220
    radius = d / 2.0

    # wind speed
    u = 10
    
    # TI
    # Since this is just to test speed, no need to have a model
    # for turbine-induced turbulence
    # Assume TI is constant throughout the farm
    ti = 6/100.0
    k = 0.003678 + 0.3837 * ti
    
    # create a square windfarm with num_turbines x num_turbines wind turbines
    xturbines, yturbines = create_farm(nx, ny, sx, sy, d)
    
    # create local points on a generic turbine rotor
    xlocal_q16, ylocal_q16, weights_q16 = create_rotor_points(d, "quad16")
    xlocal_r, ylocal_r, weights_r = create_rotor_points(d, "random")
    xlocal_c, ylocal_c, weights_c = create_rotor_points(d, "cross")

    # initiate wind speed
    speed =  zeros(len(xturbines)) + u
    
    index = -1
    time_q16 = 0
    time_r = 0
    time_c = 0
    time_anl = 0
    for i in range(nx):
        for j in range(ny):
            index += 1
            xt = xturbines[index]
            yt = yturbines[index]

            # loop over upstream turbines
            upindex = -1
            deficit_q16 = 0
            deficit_r = 0
            deficit_anl = 0
            deficit_c = 0
            num_ups = 0
            for iup in range(i):
                for jup in range(ny):
                    upindex += 1
                    xup = xturbines[upindex]
                    yup = yturbines[upindex]
                    if xup > xt: continue
                    dx = xt - xup
                    dy = abs(yt - yup)

                    ct = interpolate_ct(speed[upindex])
                    epsilon = 0.2 *   sqrt(0.5*((1+  sqrt(1-ct))/  sqrt(1-ct)))
                    sigma = k*dx + epsilon * d
                    c = 1-  sqrt(1- clip(ct/(8*(sigma/d)**2),0,1))

                    if dy / sigma > 5: continue
                    num_ups += 1
                    sigmasqr = sigma**2
                    
                    # All previous lines are shared by the two approaches
                    # Hence, they are not included in the time measurement
                    
                    ## Analytical approach
                    start_anl = time.process_time()
                    for t in range(epoch):
                        rho = dy/radius
                        sigma_hat = (sigma/radius)**2 # I've squared this as is the value used all the time
                        psi = calc_psi(rho, sigma_hat, 20)
                        deficit_anl += 2* speed[upindex] * c * sigma_hat *  exp(-(1+rho**2)/(2*sigma_hat)) * psi
                    deficit_anl /= epoch
                    end_anl = time.process_time()
                    time_anl += end_anl - start_anl

                    ## Numerical approach (Quad 16)
                    start_num = time.process_time()
                    for t in range(epoch):
                        points = create_turbine_points(xlocal_q16,ylocal_q16,yt,yup)
                        dr2 =  power(norm(points, axis = 1),2)
                        deficit_q16 +=  sum(weights_q16 * speed[upindex] * c *  exp(- divide(dr2, 2.0 * sigmasqr)))
                    deficit_q16 /= epoch
                    end_num = time.process_time()
                    time_q16 += end_num - start_num
                    
                    ## Numerical approach (Random 100)
                    start_num = time.process_time()
                    for t in range(epoch):
                        points = create_turbine_points(xlocal_r,ylocal_r,yt,yup)
                        dr2 =  power(norm(points, axis = 1),2)
                        deficit_r +=  sum(weights_r * speed[upindex] * c *  exp(- divide(dr2, 2.0 * sigmasqr)))
                    deficit_r /= epoch
                    end_num = time.process_time()
                    time_r += end_num - start_num
                    
                    ## Numerical approach (cross)
                    start_num = time.process_time()
                    for t in range(epoch):
                        points = create_turbine_points(xlocal_c,ylocal_c,yt,yup)
                        dr2 =  power(norm(points, axis = 1),2)
                        deficit_c +=  sum(weights_c * speed[upindex] * c *  exp(-  divide(dr2, 2.0 * sigmasqr)))
                    deficit_c /= epoch
                    end_num = time.process_time()
                    time_c += end_num - start_num

            if num_ups > 0: speed[index] = (u - deficit_anl / num_ups)

    print(str(i) + " " + float_to_string(deficit_anl) + " " + float_to_string(deficit_q16) + " " +
          float_to_string(deficit_c) + " " + float_to_string(deficit_r))
    #print(str(nx*ny) + " " + float_to_string(time_anl) + " " + float_to_string(time_q16) + " "
    #       + float_to_string(time_r))
    print(str(nx*ny) + " " + float_to_string((1-time_anl/time_r)*100) + " " + 
          float_to_string((1-time_q16/time_r)*100))
    print("TURBINES     TIME (ANL)    TIME (NUM)    DIFF PERC")
    #print(str(num_turbines) + "                " + str(time_anl) + "       " + str(time_num) + "       " + str((time_anl-time_num)/time_anl*100))

def float_to_string(f):
    return str(round(f,3))

def calc_psi(rho,sigsqr,ns): # sigma here is now sigma**2
    i00 = i0(rho/sigsqr)
    i11 = i1(rho/sigsqr)
    fkm1 = 1
    gkm1 = 0
    alpha = rho**2 / sigsqr
    sum0 = 0
    sum1 = 0
    for kk in range(ns):
        k = kk + 1
        fk = (fkm1 + alpha*gkm1)/k
        gk = (fk + 2*gkm1)/(2*k)
        rr = (0.5/sigsqr) ** k
        new0 = fk * rr
        new1 = gk * rr
        sum0 += new0
        sum1 += new1
        if (abs(new0 / sum0) <= 0.001) and (abs(new1 / sum1) <= 0.001): break
        fkm1 = fk
        gkm1 = gk
    return i00 * sum0 - rho/sigsqr * i11 * sum1

def generic_Phi(x, i00, ki11, tau, ns):
    fkm1 = 1
    gkm1 = 0
    sum0 = 0
    sum1 = 0
    for kk in range(ns):
        k = kk + 1
        fk = (fkm1 + tau*gkm1)/k
        gk = (fk + 2*gkm1)/(2*k)
        rr = x ** k
        new0 = fk * rr
        new1 = gk * rr
        sum0 += new0
        sum1 += new1
        if (abs(new0 / sum0) <= 0.005) and (abs(new1 / sum1) <= 0.005): break
        fkm1 = fk
        gkm1 = gk
    return i00 * sum0 - ki11 * sum1
                    
def interpolate_ct(u):
    cts = [0.0000,0.8175,0.7921,0.7864,0.7888,0.7907,0.7920,0.7918,0.7903,0.7882,0.7858,0.7833,
        0.7785,0.7785,0.7785,0.7785,0.7785,0.7785,0.7785,0.7785,0.7785,0.7785,0.7785,0.7785,0.7785,
        0.7785,0.7785,0.7785,0.7785,0.7815,0.7589,0.6144,0.4986,0.4163,0.3519,0.2998,0.2569,0.2213,
        0.1915,0.1664,0.1452,0.1273,0.1120,0.0990,0.0879,0.0783,0.0701,0.0629,0.0568,0.0514,0.0468,
        0.0,0.0]

    us = [0.00,3,3.54,4.06,4.55,5.00,5.42,5.80,6.15,6.46,6.73,6.96,7.15,7.31,7.42,7.50,7.53,7.54,
        7.58,7.67,7.80,7.97,8.17,8.42,8.70,9.02,9.38,9.78,10.2,10.6,10.8,11.1,11.6,12.2,12.8,13.4,14.1,
        14.7,15.4,16.1,16.9,17.6,18.4,19.2,20.0,20.8,21.6,22.4,23.3,24.1,25,25.02,50.0]

    if u < us[0]:
        return cts[0]
    elif u > us[-1]:
        return cts[-1]
    else:
        f = interp1d(us, cts)
        return f(u).max()
    
def uniform_weight(np):
    return [1.0 / np for k in range(np)]
    
def create_rotor_points(d, flag, npoints = 0):
    if flag == "quad16":
        rs = [0.444,0.230,0.444,0.230,0.444,0.230,0.444,0.230,0.444,
                0.230,0.444,0.230,0.444,0.230,0.444,0.230]
        ts = [0,0.393,0.785,1.178,1.571,1.963,2.356,2.749,3.142,3.534,
                3.927,4.320,4.712,5.105,5.498,5.890]
        x = [d*r*  cos(t) for r, t in zip(rs,ts)]
        y = [d*r*  sin(t) for r, t in zip(rs,ts)]  
        w = uniform_weight(len(rs))
        
    elif flag  == "random":
        x = []
        y = []
        for k in range(100):
            r =   sqrt(random.uniform(0,1)) * d / 2.0
            th = random.uniform(0,1) * 2 *   pi
            x.append(r *   cos(th))
            y.append(r *   sin(th))
        w = uniform_weight(len(x))

    elif flag == "cross":
        xs = [0.047,0.142,0.237,0.332,0,0,0,0,-0.048,-0.144,-0.239,
                -0.334,0,0,0,0]
        ys = [0,0,0,0,0.048,0.143,0.238,0.334,0,0,0,0,-0.048,
                -0.143,-0.239,-0.334]
        x = [xss * d for xss in xs]
        y = [yss * d for yss in ys]
        w = uniform_weight(len(xs))
        
    elif flag == "flower":
        b = round(2 *  sqrt(npoints))  # number of boundary points
        phi = ( sqrt(5) + 1) / 2  # golden ratio
        x = []
        y = []
        for k in range(1, npoints + 1):
            r = radius(k, npoints, b) * d
            theta = 2 *  pi * k / phi**2
            x.append(r *  cos(theta) / 2)
            y.append(r *  sin(theta) / 2)
        w = uniform_weight(len(x))
        
    return  array(x),  array(y),  array(w)

def radius(k, n, b):
    if k > n - b:
        return 1  # put on the boundary
    else:
        return sqrt(k - 1 / 2) / sqrt(n - (b + 1) / 2)  # apply square root
 
def create_farm(nx, ny, sx, sy, d):
    xs =  zeros(nx*ny)
    ys =  zeros(nx*ny)
    k = -1
    for i in range(nx):
        x = i * sx * d
        for j in range(ny):
            k += 1
            y = j * sy * d
            xs[k] = x
            ys[k] = y     
    return xs, ys

################################################################################

#test_speed_nonsym(nx = 25, ny = 25, sx = 5, sy = 3, epoch = 5, nps_in = 16, flag = 1)
#test_speed_sym(num_turbines = 30, sx = 5, sy = 3, epoch = 1)
    
#for case in ["linear", "rms", "lm"]:
#    for yy in [0.0001, 1]:
#        for anl in ["circ", "rect"]:
#            solve_HornsRev("/home/hydro/karim/AWM/HornsRev", "windturbines.txt",
#                        "wind-turbine-1.tbl", "yaw_angles_horns_rev_270.csv",
#                        case, yy, anl)

#create_figure("/home/hydro/karim/AWM/HornsRev", "wind_speeds_linear",
#              "wind_speeds_rms", "wind_speeds_lm" ,
#              "yaw_angles_horns_rev_270.csv", "windturbines.txt", 
#              "wind-turbine-1.tbl")

#create_figure_no_yaw("/home/hydro/karim/AWM/HornsRev", "wind_speeds_linear",
 #             "wind_speeds_rms", "wind_speeds_lm" ,
 #             "yaw_angles_horns_rev_270.csv", "windturbines.txt", 
 #             "wind-turbine-1.tbl")

#create_figure_no_yaw_rect("/home/hydro/karim/AWM/HornsRev",
#                "yaw_angles_horns_rev_270.csv", "windturbines.txt", 
#                "wind-turbine-1.tbl")

