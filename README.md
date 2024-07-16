# Rotor-averaged wind-speed deficit of a turbine in a non-axisymmetric Gaussian wake

This repository contains a Python implementation of an analytical expression for the rotor-averaged wind-speed deficit of a turbine placed in a non-axisymmetric Gaussian wake. 

## Background

![wake](./wake.png)

For a set of axes $y_n$-$z_n$ placed at the center of the wake, the normalised wind-speed deficit here is defined as
$$W = C\,e^{-(y_n+\omega z_n)^2/(2\sigma_y^2)} e^{-z_n^2/(2\sigma_z^2)}$$
where, 
- $C$ is a streamwise scaling function
- $\omega$ is a wind-veer coefficient that relates to the difference in wind direction ($\Delta \alpha$) across the top and bottom tips of a the wake source via $\omega=\Delta \alpha (x/D)$ with $x$ being the distance downstream of the wake source and $D$ being the wake-source diameter
- $\sigma_y$ and $\sigma_z$ are the wake standard deviations in the $y_n$ and $z_n$ directions, respectively.

The presented expression integrates the equation for $W$ across a circular disk, depicting the rotor of a downstream turbine, which is offset from the wake center by the radial distance $\rho$ and the angle $\delta$. Since, $\sigma_y$ and $\sigma_z$ are not equal in the general case that the wake source is a yawed turbine, the two standard deviations are related via the eccentricity $\xi$ such that $$\xi^2=1-(\sigma_y/\sigma_z)^2$$

Here, we use $\sigma_z=\sigma$ and hence $\sigma_y = \sigma \sqrt{1-\xi^2}$. The rotor-averaged value of $W$ is
$$\tilde{W} \approx 2C e^{-\rho^2/(2\hat{\sigma}^2)}\left(\mu_0 \left(1+2\mathcal{P}_{\textrm{ns}} \right)-\frac{2\sigma_*^2}{R^2} e^{-R^2/(2\sigma_*^2)} \mathcal{P}_{\textrm{ns}} \left[\frac{\lambda}{\rho} I_{1}\left(\frac{R\rho}{\sigma_\textrm{s}^2}\right)+\frac{\lambda^2}{\rho^2} I_{2}\left(\frac{R\rho}{\sigma_\textrm{s}^2}\right)\right]\right)$$

where $R$ is the radius of the turbine whose rotor-averaged deficit is sought. Multiple constants are included in this expression and are defined as
- $\tan{\phi_\textrm{ns}} = 2\omega / (\xi^2 - \omega^2)$
- $\tan{\phi_\textrm{s}} = \omega + \tan{\delta}\left(1-\xi^2\right)  / (1+\omega\tan{\delta})$
- $\phi=2\phi_\textrm{s}-\phi_\textrm{ns}$

- $\sigma_\textrm{s}^2 = 
    \sigma^2\left(1-\xi^2\right) \cos{\phi_\textrm{s}}\, / \, (\cos{\delta}+\omega\sin{\delta})$
- $\sigma_\textrm{ns}^2
    =2\sigma^2(1-\xi^2) / \sqrt{(\omega^2-\xi^2)^2+4\omega^2}$
- $\sigma_*^2
    =2\sigma^2(1-\xi^2) / (2+\omega^2-\xi^2)$
- $\hat{\sigma}^{-2} = \sigma_*^{-2} + \sigma_\textrm{ns}^{-2}\cos{(2\delta-\phi_\textrm{ns})}$
- $\lambda=R\sigma_\textrm{s}^2/\sigma_*^2$
- $\mathcal{P}_{\textrm{ns}} = e^{-\chi_{\textrm{ns}}^2 \cos{\phi}} \cos{(\chi_{\textrm{ns}}^2 \sin{\phi})}-1$
- $\chi_{\textrm{ns}} = \rho \sigma_*^2/(2\sigma_\textrm{ns} \sigma_\textrm{s}^2)$

The integral $\mu_0$ is 
$$\mu_0 = \int\limits_{0}^{1} \eta e^{-\eta^2 R^2/(2\sigma_*^2)} I_0 \left(\frac{\eta R\rho}{\sigma_\textrm{s}^2}\right)~d\eta$$
and is evaluated as
$$\mu_0 = \frac{\sigma_{*}^2}{R^2} e^{-R^2/(2\sigma_{*}^2)} \Psi(R,\rho,\sigma_\textrm{s},\sigma_*)$$
such that
$$\Psi(R,\rho,\sigma_\textrm{s},\sigma_*) = 
    I_0\left(\frac{R\rho}{\sigma_\textrm{s}^2}\right) \sum_{k\ge1}  \left[\left(\frac{R^2}{2 \sigma_{*}^2}\right)^{k} f_k(\tau^2) \right] - \frac{R\rho}{\sigma_\textrm{s}^2}\, I_1 \left(\frac{R\rho}{\sigma_\textrm{s}^2}\right) \sum_{k\ge1} \left[\left(\frac{R^2}{2 \sigma_{*}^2}\right)^{k} g_k(\tau^2)\right]$$
where $\tau = \rho\sigma_*/\sigma_\textrm{s}^2$ and $$f_k(\tau) = \frac{f_{k-1}(\tau) + \tau g_{k-1}(\tau)}{k}, \quad g_k(\tau) = \frac{f_{k}(\tau) + 2g_{k-1}(\tau)}{2k}$$
with the initial conditions $f_0=1$, $g_0=0$.

## Usage

To use this function, include the following lines in your code:
```
from REWS import anl_rews
w_c = anl_rews(R, rho, delta, xi, omega)
```

The function "anl_rews" takes in the radius "R" and the radial offset "rho" in a non-dimensional form by normalising the dimensional quantities by the wake standard deviation sigma. It then returns the rotor-averaged normalised deficit $\tilde{W}$ referenced to the streamwise scaling function $C$ (i.e., $\tilde{W}/C$).
