import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, solve_ivp

h = 0.67 

Om_m = 0.315            # matter density
Om_b = 0.022/(h*h)      # baryonic matter density
Om_lambda = 0.685       # vacuum energy density

rho_c = 2.78e11 *h*h    # in units of M_s/Mpc^3

# mean matter density of the universe
rho_0 = Om_m * rho_c    # in units of M_s/Mpc^3
delta_c = 1.686         # Critical density for spherical collapse

# for normalising the power spectrum
sigma8_obs = 0.811      # Observed value of sigma_8
n = 0.965               # power spectrum P(k)~ k^n

z_range = np.linspace(0.1,20, 100)

from numba import jit

f_b = Om_b/Om_m     # fraction of baryonic matter
H0 = h/9.78e9         # in units of 1/yr

@jit(nopython=True)
def H(z): #in km/yr*km
    return H0 * np.sqrt( Om_m*(1+z)**3 + Om_lambda )

def W(k,R):
    kR = k*R
    return 3*(np.sin(kR) - kR*np.cos(kR))/(kR**3)

def P(k,z,cutoff=np.inf):
    if k<=cutoff:    
        return Pi(k)*(T(k)*D(z))**2
    else:
        return 0

# Primordial power spectrum (un-normalised)
def Pi(k):
    return k**n

# Transfer Function
def T(k):
    Gamma = Om_m * h * np.exp( -Om_b* (1+np.sqrt(2*h)/Om_m) )
    q = k / (Gamma*h)
    value = 1 + (3.89*q) + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4
    return np.log(1+2.34*q) / (2.34*q*value**(1/4))

def D(z):
    Dz = (Om_m + 0.4545*Om_lambda) / (Om_m*(1+z)**3 + 0.4545*Om_lambda)
    return Dz**(1/3)

M_range = np.logspace(6, 16, 500)  # Mass range in solar masses

def poly_fit(x_arr,y_arr,degree,plot=0):
    # Perform polynomial fit
    coefficients = np.polyfit(x_arr, y_arr, degree)

    # Generate the fitted polynomial function
    polynomial = np.poly1d(coefficients)

    # Plot the original data and the fitted polynomial
    if (plot):
        plt.figure(figsize=(12, 7))
        plt.plot(x_arr, polynomial(x_arr), 'r', lw=1.5, label=f'Polynomial Fit (Degree {degree})')
        plt.plot(x_arr, y_arr, '--',color='k', label='Original Data')
        plt.xscale('log')
        plt.yscale('log')  
        plt.xlabel('$\log M$')
        plt.ylabel('$\sigma_o(M)$')
        plt.legend()
        plt.title('Polynomial Fit')
        plt.show()
    
    return coefficients

degree = 9

logM = np.log(M_range)

R8 = 8/h      # units of Mpc

def sigma_R(R,z=0,cutoff=1000):    
    # Integrand for the variance calculation
    def integral(k):
        return P(k,z,cutoff) * k**2 * W(k,R)**2 / (2 * np.pi**2)
    
    sigma2, _ = quad(integral, 0, np.inf, limit=500, epsabs=1e-5, epsrel=1e-5)
    return np.sqrt(sigma2)

def sigma_M(M,z=0,cutoff=1000):
    # Compute R corresponding to the mass M
    R = (3 * M / (4 * np.pi * rho_0))**(1/3)
    
    return sigma_R(R,z,cutoff)

A = sigma8_obs/sigma_R(R8,0)

sigma_range = [A*sigma_M(M) for M in M_range]

cof = poly_fit(logM,sigma_range,degree=9)
d_cof = np.polyder(cof)

fit_deriv = np.poly1d(d_cof)
    
def PS_MassFunc(M, z):
    # Generate the fitted polynomial function
    polynomial = np.poly1d(cof)
    sigma0 = polynomial(np.log(M))
    d_sigma = abs(fit_deriv(np.log(M)))
    factor = np.sqrt(2 / np.pi) * delta_c/(sigma0**2 * D(z)) * abs(fit_deriv(np.log(M)))
    exponent = (-delta_c**2 / (2 * sigma0**2 * D(z)**2))
    return factor * np.exp(exponent)

def Madau(z):
    return 0.01 * (1+z)**2.6 / ( 1+((1+z)/3.2)**6.2 )

def Harikane(z):
    return 1 / ( 61.7 * (1+z)**(-3.13) + 10**(0.22*(1+z)) + 2.4 * 10**(0.5*(1+z)-3) )

def New_SFR(z): #Khaire
    return 10 ** (-2) * (2.01 + 8.48 * z) / (1 + (z / 2.5) ** 3.09)

mean_cx = np.linspace(2.6 * 10 ** (39), 3.7 * 10 ** (39)) #erg s^-1 / M yr^-1, Dijkstra (2012)
mean_cx_scat = mean_cx * np.e ** (1 / 2 * (0.4 * np.log(10)) ** 2)

Om_m_new = 0.27
Om_lambda_new = 0.73
h_new = 0.7

def Hopkin_SFRD(z): #Hopkin and Beacom (2006)
    a = 0.017
    b = 0.13
    c = 3.3
    d = 5.3
    return (a + b * z) * h_new / (1 + (z / c) ** d) #in M yr^-1 cMpc^-3

n_gamma = 4800 #Number of ionizing photons produced per baryon (proton), dimensionless
n_gamma2 = 7780
f_esc = 0.1 #Fraction of photons escaping the star forming halo
f_esc2 = 0.2
m_p = 1.67 * 10 ** (-27) * 5.03 * 10 ** (-31) #mass of proton in Solar Masses

def RK4(f, x_0, y_0, h, x_n):
    steps = int(np.abs(np.abs(x_n - x_0) / h))
    x_arr = np.zeros(steps + 1)
    y_arr = np.zeros_like(x_arr)
    
    x_arr[0] = x_0
    y_arr[0] = y_0
    counter = 0

    for i in range(1, steps + 1):
        counter += 1
        k_1 = f(x_arr[i - 1], y_arr[i - 1])
        k_2 = f(x_arr[i - 1] + h / 2, y_arr[i - 1] + h * k_1 / 2)
        k_3 = f(x_arr[i - 1] + h / 2, y_arr[i - 1] + h * k_2 / 2)
        k_4 = f(x_arr[i - 1] + h, y_arr[i - 1] + h * k_3)

        y_arr[i] = y_arr[i - 1] + h * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        x_arr[i] = x_arr[i - 1] + h

        if y_arr[i] >= 1:
            for j in range(counter, steps + 1):
                y_arr[j] = 1
                x_arr[j] = x_arr[j - 1] + h
            break
    return x_arr, y_arr

X_H = 0.75 #Fraction of hydrogen
alpha_B = 2.59 * 10 ** (-13) #Case B recombination coefficient at T=3*10^4 K in cm^3 * s^-1

alpha_B_unit = alpha_B * (3.24 * 10 ** (-25)) ** 3 / (3.17 * 10 ** (-8)) #Coefficient in Mpc^3 * yr^-1

def n_H(z): #Mean proper number denisty of hydrogen atoms in Mpc^-3
    return X_H * Om_b * rho_c * (1 + z) ** 3 / m_p

def dt_dz(z): #In yr
    return -1 / ((1 + z) * H(z))

def dN_dt_Mad(z):
    return Madau(z) * (1 + z) ** 3 * n_gamma2 * f_esc / m_p

def Clump(z):
    #if z >= 6:
     #   return 1 + 9 * (7 / (1 + z)) ** 2
    #else:
     #   return 10
    #return 3
    return 2.9 * ((1 + z) / 6) ** (-1.1)

def RHS_Mad(z, f_HII):
    return dN_dt_Mad(z) * dt_dz(z) / n_H(z) - alpha_B_unit * n_H(z) * f_HII * Clump(z) * dt_dz(z)

Mad_data = RK4(RHS_Mad, 20, 0, -0.1, 0.01)