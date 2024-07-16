import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from REWS import anl_rews

def verify(daoa_deg, yaw_deg):
    matplotlib.use('Agg')
    plt.rcParams.update({"text.usetex": True, "font.family": "sans-serif"})

    # Create figure
    fig, axs = plt.subplots(ncols = 4, nrows = 4, figsize=(10,8))
    
    # Convert to radians
    daoa = daoa_deg * np.pi / 180.0
    yaw = yaw_deg * np.pi / 180.0

    # beta here is rho/sigma
    betas = [0.001,0.5,1,1.5,2,2.5,3,3.5,4]
    betas_refined = np.linspace(0.0001, 4, 50, endpoint=True)
    
    # different values for x/D
    xds = [4,6,8,10]
    
    # Different values for the angle delta
    deltas = [0, np.pi/4, 3*np.pi/4]

    # Figure related
    markers = ["o","^","s","v","D"]
    colors = ["k", "tab:red", "tab:blue", "tab:orange"]
    letters = ["a","b","c","d","e","f","g","h","i","j",
        "k","l","m","n","o","p","q","r","s","t","u",
        "v","w","x","y","z"]
    delta_names = ["0", "\\pi/4", "3\\pi/4"]
    handles = []
    
    # The radius of NREL-5MW
    # This value is redundant
    radius = 124.6 / 2.0 # NREL-5MW
    
    # Use 2000 points following a sunflower distrbution as a reference
    # for the analytical expression being verified
    pss = sunflower(radius, n = 2000)
    
    row_index = 0
    column_index = -1

    for ct in [0.4,0.8]:
        for ti in [0.05, 0.12]:
            for xd in xds:
                gamma, xi, omega = calc_constants(ct, ti, yaw, xd, daoa)

                column_index += 1
                if column_index > 3:
                    row_index += 1
                    column_index = 0

                ax = axs[row_index, column_index]

                for delta in deltas:
                    idelta = deltas.index(delta)
                    label = "$\\delta=" + delta_names[idelta] + "$"
                    
                    curve,  = ax.plot(betas_refined, [anl_rews(gamma, beta, delta, xi, omega) for beta in betas_refined],
                        linewidth = 1, color = colors[idelta], label= label, zorder = 1)
                    
                    ax.scatter(betas, [solve_numerical(pss, radius / gamma, beta, delta, xi, omega) for beta in betas],
                        marker=markers[idelta], s = 16, facecolors='white', edgecolor = colors[idelta],
                        zorder = 20, label= label, linewidth = 0.75)
                    
                    if row_index == column_index == 0: handles.append(curve)

                    ax.tick_params(axis='both', which='major', labelsize=8)
                    ax.text(0.22,0.88,"$x/D=" + str(xd) + "$",transform=ax.transAxes,fontsize = 10)
                    ax.text(0.63,0.78,"$\\omega=" + str(round(omega,2)) + "$",transform=ax.transAxes,fontsize = 10)
                    ax.text(0.63,0.66,"$R/\\sigma=" + str(round(gamma,1)) + "$",transform=ax.transAxes,fontsize = 10)
                    ax.text(0.63,0.9,"$\\xi=" + str(round(xi,2)) + "$",transform=ax.transAxes,fontsize = 10)
                    ax.text(0.1,0.88,"(" + letters[column_index + 4*row_index] + ")",transform=ax.transAxes,fontsize = 10)
                    ax.text(0.02,0.18,"$C_t=" + str(ct) + "$",transform=ax.transAxes,fontsize = 9)
                    ax.text(0.02,0.06,"$T_i=" + str(int(ti*100)) + "\%$",transform=ax.transAxes,fontsize = 9)
                    
                    ax.set_ylim([-0.05,1.02])
                    ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
                    if column_index == 0: ax.set_ylabel("$\\tilde{W}/C$", fontsize = 9)
                    else: ax.yaxis.set_ticklabels([])
                    
                    ax.set_xlim([-0.2,4.2])
                    ax.set_xticks([0,1,2,3,4])
                    if row_index == 3: ax.set_xlabel("$\\rho/\\sigma$", fontsize = 9)
                    else: ax.xaxis.set_ticklabels([])
                         
    fig.legend(handles = handles, loc = [0.33, 0.96], ncol=4, prop={'size': 10}, frameon=False)
    fig.subplots_adjust(bottom=0.075, top=0.955, left=0.075, right=0.975, wspace=0.1, hspace=0.1)
    fig.savefig("verify.png", dpi=300)

def sunflower(radius, n = 1000):
    alpha = 2
    b = round(alpha * np.sqrt(n))  # number of boundary points
    phi = (np.sqrt(5) + 1) / 2  # golden ratio
    rs= []
    ts = []
    for k in range(1, n + 1):
        r = calc_radius(k, n, b) / 2.0
        theta = 2 * np.pi * k / phi**2
        rs.append(r)
        ts.append(theta)
    return [2*radius*r*(np.cos(t) + 1j * np.sin(t)) for r, t in zip(rs,ts)]

def calc_radius(k, n, b):
    if k > n - b:
        return 1  # put on the boundary
    else:
        return np.sqrt(k - 1 / 2) / np.sqrt(n - (b + 1) / 2) 

def calc_constants(ct, ti, yaw, xd, daoa):
    eps = 0.2 * np.sqrt( ( 1.0 + np.sqrt(1.0 - ct) ) / ( 2 * np.sqrt(1.0 - ct) ) )
    ks = 0.003678 + 0.3837 * ti
    sz0 = np.sqrt( ( 1.0 + np.sqrt(1.0 - ct * np.cos(yaw)) ) / ( 8.0 * ( 1 + np.sqrt(1-ct) ) ) )
    gamma = 0.5 / (ks * xd + eps)
    omega = daoa * xd
    syd = ks * xd + np.cos(yaw) * sz0
    szd = ks * xd + sz0
    xi = np.sqrt(1-(syd / szd)**2)
    return gamma, xi, omega

def solve_numerical(pss, sigma, beta, delta, xi, omega):
    sigma_z = sigma
    sigma_y = sigma * np.sqrt(1-xi**2)
    turb_center = beta * sigma * np.cos(delta) + 1j * beta * sigma * np.sin(delta)
    ps = [p + turb_center for p in pss]
    ws = [np.exp(-(p.real + omega*p.imag)**2 / 2.0 / sigma_y**2) * np.exp(-(p.imag)**2 / 2.0 / (sigma_z**2)) for p in ps]
    return np.mean(ws)

verify(7, 20)