import numpy as np
from scipy.optimize import curve_fit
import argparse
import sys

def get_args():
    """Define the task arguments with the default values.                         
    Returns:                                                                      
      experiment parameters                                                     
    """
    
    args_parser = argparse.ArgumentParser()
    
    args_parser.add_argument(
        '--qmc_ks',
        help='Filename of the QMC twisted energies in HA.',
        type=str,
        default='E_qmc_ks.txt'
    )

    args_parser.add_argument(
        '--dft_ks',
        help='Filename of the DFT twisted energies in Ry.',
        type=str,
        default='E_dft_ks.txt'
    )

    args_parser.add_argument(
        '--dft_dense',
        help='Value of the dft energy with dense k-grid in Ry',
        type=float,
        default=0.0
    )

    args_parser.add_argument(
        '--atsup_dense',
        help='Difference of the atsup energies epsilon_ks - epsilon_dense in Ha.',
        type=float,
        default=0.0
    )
    
    return args_parser.parse_args()

def Eks(Edft, E_ta, b):
    return E_ta+b*Edft

def main():
    args=get_args()

    # Inputs: Energies of QMC, DFT, and DFT dense. In addition, positron energy
    #         is given as the differense E_atsup_loose-E_atsup_dense.

    Eqmc=np.loadtxt(args.qmc_ks)
    Edft=0.5*np.loadtxt(args.dft_ks)
    
    if(Eqmc.shape[0]!=Edft.shape[0]):
        sys.exit("QMC and DFT have different number of twists.")

    Edft_dense=0.5*args.dft_dense
    Epositron=args.atsup_dense

    ydata=Eqmc[:,0]
    xdata=Edft-Edft_dense+Epositron

    popt,pcov=curve_fit(Eks,xdata,ydata)

    ha2ev=27.2114
    
    print(ha2ev*popt)
    print("-")
    print(ha2ev*pcov)

if __name__=='__main__':
    main()
