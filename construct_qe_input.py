from ase.build import bulk
from ase.build import find_optimal_cell_shape, get_deviation_from_optimal_cell_shape
from ase.build import make_supercell
from ase.visualize import view
from ase.calculators.espresso import Espresso
from ase.constraints import UnitCellFilter
from ase.optimize import LBFGS
from ase.io.espresso import write_espresso_in
import sys
import numpy as np
import argparse

def get_args():
  """Define the task arguments with the default values.                         
  Returns:                                                                      
      experiment parameters                                                     
  """

  args_parser = argparse.ArgumentParser()

  # Data files arguments                                                        
  args_parser.add_argument(
    '--atom',
    help='Chemical element of an atom. Currently implemented: Si, C.',
    type=str,
    required=True
  )

  args_parser.add_argument(
    '--lattice_constant',
    help='Lattice constant in Bohr.',
    type=float,
    required=True
  )

  args_parser.add_argument(
    '--structure',
    help='''
    The crystal structure: either sc, fcc, bcc, hcp, 
    diamond, zincblende, rocksalt, cesiumchloride, fluorite or wurtzite.
    ''',
    type=str,
    default='diamond'
  )

  args_parser.add_argument(
    '--Natoms',
    help='Number of atoms in the supercell',
    type=int,
    default=0
  )

  args_parser.add_argument(
    '--vacancy',
    help='Number of atoms to pop from the beginning of atom list.',
    type=int,
    default=0
  )

  args_parser.add_argument(
    '--remove_electrons',
    help='Number of electrons to pop from the calculations.',
    type=int,
    default=0
  )

  args_parser.add_argument(
    '--ecutwfc',
    help='Kinetic energy cutoff (Ry) for wavefunctions',
    type=float,
    required=True
  )

  args_parser.add_argument(
    '--ecutrho',
    help='Default: 4*ecutwfc',
    type=float,
    default=-1.0
  )

  args_parser.add_argument(
    '--input_dft',
    help='The correlation functional',
    type=str,
    default='PBE'
  )

  args_parser.add_argument(
    '--kpts',
    help='The number of kpoints at each direction, separated by comma and NO spaces, like: 2,2,2',
    type=str,
    default='1,1,1'
  )

  args_parser.add_argument(
    '--conv_thr',
    help='Convergence threshold for selfconsistency.',
    type=float,
    default=1e-9
  )

  args_parser.add_argument(
    '--mixing_beta',
    help='The mixing parameter.',
    type=float,
    default=0.7
  )

  args_parser.add_argument(
    '--view',
    help='Visualize the structure.',
    type=int,
    default=0
  )
  
  return args_parser.parse_args()

def construct_supercell(conf, N):

    P1 = find_optimal_cell_shape(conf.cell, N, 'sc')
    sc = make_supercell(conf,P1)

    return sc

def setup_dft_input(args):

  # Define density cutoff
  if args.ecutrho < 0.0:
    ecutrho  = 4 * args.ecutwfc
  else:
    ecutrho=args.ecutrho

  # prefix
  name=args.atom
  if(name=='Si'):
    prefix='silicon'
  elif(name=='C'):
    prefix='diamond'
  else:
    sys.exit('Unknown atoms provided.')

  # The system charge
  tot_charge=args.remove_electrons

  input_data = {
    'control':{
      'prefix': prefix,
      'restart_mode': 'from_scratch',
      'calculation': 'scf',
      'wf_collect': True,
      'pseudo_dir':'./',
      'outdir': './'},
    'system': {
      'ibrav': 0,
      'ecutwfc': args.ecutwfc,
      'ecutrho': ecutrho,
      'input_dft': args.input_dft,
      'nosym': True,
      'noinv': True,
      'tot_charge': tot_charge},
    'electrons':{
      'diagonalization':'cg',
      'conv_thr': args.conv_thr,
      'mixing_mode': 'plain',
      'mixing_beta': args.mixing_beta}
  }
  
  return input_data
  
def main():

    args = get_args()

    # Define input file parameters
    input_data = setup_dft_input(args)

    # Define k point grid
    kpts_string=args.kpts.split(',')
    kpts=(int(kpts_string[0]),int(kpts_string[1]),int(kpts_string[2]))
    
    pseudopotentials = {'C': 'c_df.UPF',
                        'Si': 'si_df.UPF'}

    conf = bulk(args.atom,
                crystalstructure=args.structure,
                a=args.lattice_constant)
    # If the required number of atoms is specified, construct a supercell
    if(args.Natoms>0):
      conf = construct_supercell(conf,args.Natoms)
      for i in range(args.vacancy):
        conf.pop()
      
    with open('in.pwscf','w') as f:
      write_espresso_in(f,conf,
                        input_data=input_data,
                        pseudopotentials=pseudopotentials,
                        kpts=kpts)

    print(conf.get_positions())

    if(args.view==1):
      view(conf)
    
if __name__ == '__main__':
    main()
