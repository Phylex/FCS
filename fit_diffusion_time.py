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
    r0 = radius_of_volume(popt[0], D)
    dpr = np.abs(r0 - radius_of_volume(popt[0] + np.sqrt(pcov[0][0]), D))
    dnr = np.abs(r0 - radius_of_volume(popt[0] - np.sqrt(pcov[0][0]), D))
    V = Veff(r0, params['z_r_ratio']*r0)
    dpv = np.abs(V - Veff(r0+dpr, (r0+dpr)*params['z_r_ratio']))
    dnv = np.abs(V - Veff(r0-dnr, (r0-dnr)*params['z_r_ratio']))
    C = Concentration(V, params['Nin'])
    dpc = np.abs(C - Concentration(V+dpv, params['Nin']))
    dnc = np.abs(C - Concentration(V-dnv, params['Nin']))
    print("Calculated values from fit result:")
    print("r_0 = {} (+{}/-{}) m".format(r0, dpr, dnr))
    print("Veff = {} (+{}/-{}) m^3".format(V, dpv, dnv))
    print("C = {} (+{}/-{}) 1/m^3".format(C, dpc, dnc))
    print("A plot of this fit was generatet at {}".format(figname))
    print()

def motility_function(tau, tau_d, z_r_ratio):
    return (1/(1+(tau/tau_d)))*\
           (1/np.sqrt(1+((1/z_r_ratio)**2*(tau/tau_d))))

def two_species_autocorrelation(tau, tau_d1, tau_d2, f_1, Np):
    f_2 = 1-f_1
    part_1 = f_1 * motility_function(tau, tau_d1, 5)
    part_2 = f_2 * motility_function(tau, tau_d2, 5)
    return (1/Np) * (part_1 + part_2)


def two_species_fit(data, name):
    print("Performing fit for autocorrelation function of two different species in the same sample")
    print("---------------------------------------------------------------------------------------")
    popt, pcov = curve_fit(two_species_autocorrelation, data['t'], data['c'],
                           p0=(500e-6, 10000e-6, 0.5, 0.2))
    print("Results of the fit:")
    print("tau_d1: {} +/- {} s".format(popt[0], np.sqrt(pcov[0][0])))
    print("tau_d2: {} +/- {} s".format(popt[1], np.sqrt(pcov[1][1])))
    print("f_1: {} +/- {}".format(popt[2], np.sqrt(pcov[2][2])))
    print("f_2: {} +/- {}".format(1-popt[2], np.sqrt(pcov[2][2])))
    print("Np: {} +/- {}".format(popt[3], np.sqrt(pcov[3][3])))
    print()

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
                       'Nin': 0.05, 'D': 400e-12, 'z_r_ratio': 5}
    alexa_546_params = {'Nin': 0.021, 'D': 390e-12, 'z_r_ratio': 5}

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
    print()
    two_species_fit(probe_1, "Probe 1")
    two_species_fit(probe_2, "Probe 2")

