import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys

def Veff(r_0, z_0):
    return np.pi**(3/2) * r_0**2 * z_0

def Concentration(V, Nin):
    return 1/(V*Nin)

def decay_time_constant(r_0, D):
    return r_0**2/(4*D)

def radius_of_volume(tau_d, D):
    return np.sqrt(tau_d*4*D)

def free_3D_correlation_function(tau, tau_d, Ninverse, D, z_r_ratio):
    r_0 = radius_of_volume(tau_d, D)
    Vef = Veff(r_0, r_0*z_r_ratio)
    C = Concentration(Vef, Ninverse)
    return (1/(Vef*C))*(1/(1+(tau/tau_d)))*\
           (1/np.sqrt(1+((1/z_r_ratio)**2*(tau/tau_d))))


def config_standard(D, Nin, z_r_ratio):
    """create a function that can be used to fit the radius of the confocal region
    without needing to retype the entire function.
    It also sets the correlation function"""
    return lambda tau, tau_d: free_3D_correlation_function(tau,
                                                           tau_d,
                                                           Nin,
                                                           D,
                                                           z_r_ratio)


def triplet_blinking(tau, T, tau_t):
    return (1-T)+T*np.exp(-tau/tau_t)


def combined_corellation_funciton(tau, tau_d, tau_t, Nin, T, D, z_r_ratio):
    return triplet_blinking(tau, T, tau_t) *\
           free_3D_correlation_function(tau, tau_d, Nin, D, z_r_ratio)


def config_combined(T, D, tau_t, Nin, z_r_ratio):
    return lambda tau, tau_d: combined_corellation_funciton(tau, tau_d, tau_t, Nin, T, D, z_r_ratio)

def perform_r_0_fit(name, data, func_generator, params):
    print("Perfoming r_0_fit on {}".format(name))
    print("-----------------------------")
    print("Fixed parameters of fitfunction:")
    for k in params.keys():
        print("{}: {}".format(k, params[k]))
    fit_func = func_generator(**params)
    popt, pcov = curve_fit(fit_func, data['t'], data['c'], p0=(350e-9))
    plt.semilogx(data['t'], data['c'], marker='x', linestyle='',
                 label='{}'.format(name))
    plt.plot(data['t'], fit_func(data['t'], *popt), color='blue',
             label='Fit of {}'.format(name))
    plt.grid()
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Autocorrelation c")
    figname = "r_0_fit_on_{}.pdf".format(name.replace(" ", '_'))
    plt.savefig(figname)
    plt.close()
    D = params['D']
    print("Parameters of fit:\ntau_d = {} +/- {} s".format(popt[0], np.sqrt(pcov[0][0])))
    V = Veff(popt[0], 5*popt[0])
    dpv = np.abs(V - Veff(popt[0]+np.sqrt(pcov[0][0]), 5*popt[0]+np.sqrt(pcov[0][0])))
    dnv = np.abs(V - Veff(popt[0]-np.sqrt(pcov[0][0]), 5*popt[0]-np.sqrt(pcov[0][0])))
    C = Concentration(V, D)
    dpc = np.abs(C - Concentration(V+dpv, D))
    dnc = np.abs(C - Concentration(V-dnv, D))
    r0 = radius_of_volume(popt[0], D)
    dpr = np.abs(r0 - radius_of_volume(popt[0] + np.sqrt(pcov[0][0]), D))
    dnr = np.abs(r0 - radius_of_volume(popt[0] - np.sqrt(pcov[0][0]), D))
    print("Calculated values from fit result:")
    print("r_0 = {} (+{}/-{}) m".format(r0, dpr, dnr))
    print("Veff = {} (+{}/-{}) m^3".format(V, dpv, dnv))
    print("C = {} (+{}/-{}) 1/m^3".format(C, dpc, dnc))
    print("A plot of this fit was generatet at {}".format(figname))
    print()

def two_species_autocorrelation():
    pass

def two_species_fit(data, params):
    pass

if __name__ == "__main__":
    # part 1 read in all the data:
    atto_488 = pd.read_excel("Atto488.xlsx")
    alexa_546 = pd.read_excel("Alexa546.xlsx")
    probe_1 = pd.read_excel("Probe1.xlsx")
    probe_2 = pd.read_excel("Probe2.xlsx")

    # prepare the data
    atto_488 = atto_488.rename(columns={"Corr.Time[μs]": 't', "Correlation": 'c'})
    atto_488['t'] = atto_488['t'] * 10**-6
    print(atto_488)
    alexa_546 = alexa_546.rename(columns={"Corr.Time[μs]": 't', "Correlation": 'c'})
    alexa_546['t'] = alexa_546['t'] * 10**-6
    print(alexa_546)
    probe_1 = probe_1.rename(columns={"Corr.Time[μs]": 't', "Correlation": 'c'})
    probe_1['t'] = probe_1['t'] * 10**-6
    probe_2 = probe_2.rename(columns={"Corr.Time[μs]": 't', "Correlation": 'c'})
    probe_2['t'] = probe_2['t'] * 10**-6

    # hand-write the parameters from the results file into python
    atto_488_params = {'T': 0.44, 'tau_t': 0.8*10**-6,
                       'Nin': 0.05, 'D': 400e-6, 'z_r_ratio': 5}
    alexa_546_params = {'Nin': 0.021, 'D': 390e-6, 'z_r_ratio': 5}

    # perform fits on the two datasets for part one
    print("Task 1")
    print("======")
    perform_r_0_fit("Atto 488", atto_488, config_combined, atto_488_params)
    perform_r_0_fit("Alexa 546", alexa_546, config_standard, alexa_546_params)
    print()
    print()
    # now we come to task 2
    print("Task 2")
    print("======")

