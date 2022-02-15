
import numpy as np
from classy import Class

def Get_LCDM_Pks(h, omega_b, omega_cdm, tau_reio, A_s, n_s, P_k_max, k_per_decade):
    LCDM = Class()
    
    # pass input parameters
    LCDM.set({'h':h,'omega_b':omega_b,'omega_cdm':omega_cdm, 'tau_reio':tau_reio,'A_s':A_s,'n_s':n_s})
    LCDM.set({'output':'mPk','lensing':'no','P_k_max_1/Mpc':P_k_max})
    
    # run class
    LCDM.compute()
    
    # retrieve pk's at z=0
    h = LCDM.h()
    k_vec = np.logspace(-1,np.log10(P_k_max),300) # units of h/Mpc
    LCDM_Pk_vec = np.zeros(len(k_vec))  # units of (Mpc/h)**3
    for k in range(len(k_vec)):
        LCDM_Pk_vec[k]=LCDM.pk(k_vec[k]*h,0.) * h**3
    
    # output
    return k_vec, LCDM_Pk_vec


def Get_DMEFF_Pks(h, omega_b, omega_cdm, tau_reio, A_s, n_s, P_k_max, omega_dmeff, sigma_dmeff, \
                  m_dmeff, npow_dmeff, dmeff_target, k_per_decade):
    DMEFF = Class()
    
    # pass input parameters
    DMEFF.set({'h':h,'omega_b':omega_b,'omega_cdm':1e-20, 'tau_reio':tau_reio,'A_s':A_s,'n_s':n_s})
    DMEFF.set({'omega_dmeff':omega_dmeff,'sigma_dmeff':sigma_dmeff,'m_dmeff':m_dmeff,\
               'npow_dmeff':npow_dmeff,'dmeff_target':dmeff_target, 'k_scalar_k_per_decade_for_pk':k_per_decade})
    DMEFF.set({'output':'mPk','lensing':'no','P_k_max_1/Mpc':P_k_max})
       
    # run class
    DMEFF.compute()
    
    # retrieve pk's at z=0
    h = DMEFF.h()
    k_vec = np.logspace(-1,np.log10(P_k_max),300) # units of h/Mpc
    DMEFF_Pk_vec = np.zeros(len(k_vec))  # units of (Mpc/h)**3
    for k in range(len(k_vec)):
        DMEFF_Pk_vec[k]=DMEFF.pk(k_vec[k]*h,0.) * h**3
    
    # output
    return k_vec, DMEFF_Pk_vec


def Get_WDM_Tx(m_dmeff, omega_dmeff, k_vec, h):
    
    alpha =0.049*pow(omega_dmeff/h/h/0.25, 0.11) * pow(h/0.7,1.22) * pow(1.0 / m_dmeff, 1.11)
    Tf = pow(1 + pow(alpha * k_vec[:], 2 * 1.12), -5.0 / 1.12)
    WDM_Tx_vec = Tf * Tf
    
    # output
    return k_vec, WDM_Tx_vec

def Calc_K_Max(omega_b, omega_cdm, h):
    
    # define parameters
    M_sol = 1.989e30        # kg
    M_min = 5e8 * M_sol
    G = 6.6743e-11          # N.m^2/kg^2
    H = 67.27 / 3.086e19    # 1/s
    Omega_m = (omega_b + omega_cdm)/h/h
    
    # calculate k_max
    k_max = ( np.pi * 3.086e22 * (((H**2*Omega_m)/(2*G*M_min)) ** (1/3)) ) /h  # h/Mpc}
    
    #output
    return k_max
