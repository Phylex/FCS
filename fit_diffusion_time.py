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

def free_3D_correlation_function(tau, Ninverse, D, r_0, z_0):
    Vef = Veff(r_0, z_0)
    C = Concentration(Vef, Ninverse)
    tau_d = decay_time_constant(r_0, D)
    return (1/(Vef*C))*(1/(1+(tau/tau_d)))*\
           (1/np.sqrt(1+((r_0/z_0)**2*(tau/tau_d))))


def config_standard(D, Nin):
    """create a function that can be used to fit the radius of the confocal region
    without needing to retype the entire function.
    It also sets the correlation function"""
    return lambda tau, r_0: free_3D_correlation_function(tau, Nin, D, r_0, 5*r_0)


def triplet_blinking(tau, T, tau_t):
    return (1-T)+T*np.exp(-tau/tau_t)


def combined_corellation_funciton(tau, Nin, T, tau_t, D, r_0, z_0):
    return triplet_blinking(tau, T, tau_t) *\
           free_3D_correlation_function(tau, Nin, D, r_0, z_0)


def config_combined(T, D, tau_t, Nin):
    return lambda tau, r_0: combined_corellation_funciton(tau, Nin, T, tau_t, D, r_0, 5*r_0)

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
    print("Parameters of fit:\nr_0 = {} +/- {} m".format(popt[0], np.sqrt(pcov[0][0])))
    V = Veff(popt[0], 5*popt[0])
    dpv = np.abs(V - Veff(popt[0]+np.sqrt(pcov[0][0]), 5*popt[0]+np.sqrt(pcov[0][0])))
    dnv = np.abs(V - Veff(popt[0]-np.sqrt(pcov[0][0]), 5*popt[0]-np.sqrt(pcov[0][0])))
    D = params['D']
    C = Concentration(popt[0], D)
    dpc = np.abs(C - Concentration(popt[0]+np.sqrt(pcov[0][0]), D))
    dnc = np.abs(C - Concentration(popt[0]-np.sqrt(pcov[0][0]), D))
    print("Calculated values from fit result:")
    print("Veff = {} (+{}/-{}) m^3".format(V, dpv, dnv))
    print("C = {} (+{}/-{}) 1/m^3".format(C, dpc, dnc))
    print("A plot of this fit was generatet at {}".format(figname))
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
    alexa_546 = alexa_546.rename(columns={"Corr.Time[μs]": 't', "Correlation": 'c'})
    alexa_546['t'] = alexa_546['t'] * 10**-6
    probe_1 = probe_1.rename(columns={"Corr.Time[μs]": 't', "Correlation": 'c'})
    probe_1['t'] = probe_1['t'] * 10**-6
    probe_2 = probe_2.rename(columns={"Corr.Time[μs]": 't', "Correlation": 'c'})
    probe_2['t'] = probe_2['t'] * 10**-6

    # hand-write the parameters from the results file into python
    atto_488_params = {'T': 0.44, 'tau_t': 0.8*10**-6,
                       'Nin': 0.05, 'D': 400e-6**2}
    alexa_546_params = {'Nin': 0.021, 'D': 390e-6**2}

    # perform fits on the two datasets for part one
    perform_r_0_fit("Atto 488", atto_488, config_combined, atto_488_params)
    perform_r_0_fit("Alexa 546", alexa_546, config_standard, alexa_546_params)
