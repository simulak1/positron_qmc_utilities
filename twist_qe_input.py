import numpy as np
import argparse
import sys

def get_args():
  """Define the task arguments with the default values.                         
  Returns:                                                                      
      experiment parameters                                                     
  """

  args_parser = argparse.ArgumentParser()

  # Data files arguments                                                        
  args_parser.add_argument(
      '--file',
      help='QE input file.',
      type=str,
      default='in.pwscf'
  )

  args_parser.add_argument(
    '--kpoints',
    help='The number of kpoints in each of the 3 directions',
    type=int,
    nargs=3,
    default=[2,2,2]
    )

  args_parser.add_argument(
    '--twists',
    help='The number of twists in each of the 3 directions',
    type=int,
    nargs=3,
    default=[2,2,2]
    )

  args_parser.add_argument(
      '--symfile',
      help='File containing k-point symmetries of a grid. Use with caution',
      type=str,
      default='nofile'
  )
  
  return args_parser.parse_args()

def fortran_nint(r):
  # Corresponds to fortran's nint-function.
  if(r-round(r)==0.5):
    return round(r)+1
  elif(r-round(r)==-0.5):
    return round(r)-1
  else:
    return round(r)


def parse_input(filename):
    with open(filename) as f:
        lines=f.readlines()
        for i in range(len(lines)):
            if lines[i].strip()[:9]=='celldm(1)':
              lattice_constant = lines[i].split('=')[1]
            elif(lines[i].strip()[:15]=='CELL_PARAMETERS'):
              pa=[]
              for j in range(1,4):
                words=lines[i+j].split()
                pa.append([float(words[0]),float(words[1]),float(words[2])])
    return lattice_constant, np.array(pa)
            
def reciprocal(La):
  Lb=np.zeros((3,3))
  Lb[0,:]=np.cross(La[1,:],La[2,:])/(np.dot(La[0,:],np.cross(La[1,:],La[2,:])))
  Lb[1,:]=np.cross(La[2,:],La[0,:])/(np.dot(La[1,:],np.cross(La[2,:],La[0,:])))
  Lb[2,:]=np.cross(La[0,:],La[1,:])/(np.dot(La[2,:],np.cross(La[0,:],La[1,:])))
  return Lb

def generate_kpoint_grid(nk):
  Nk=nk[0]*nk[1]*nk[2]
  kvec=np.zeros((Nk,3))
  for i in range(nk[0]):
    for j in range(nk[1]):
      for k in range(nk[2]):

       	   n = k + j*nk[2] + i*nk[1]*nk[2]
           # The zeros below are the offsets from the origin!
           kvec[n,0] = float(i)/nk[0] + 0/2/nk[0]
           kvec[n,1] = float(j)/nk[1] + 0/2/nk[1]
           kvec[n,2] = float(k)/nk[2] + 0/2/nk[2]

           kvec[n,0] = kvec[n,0]-fortran_nint(kvec[n,0])
           kvec[n,1] = kvec[n,1]-fortran_nint(kvec[n,1])
           kvec[n,2] = kvec[n,2]-fortran_nint(kvec[n,2])

           if(abs(kvec[n,0])<0.00000000000001):
             kvec[n,0]=0.0
           if(abs(kvec[n,1])<0.00000000000001):
             kvec[n,1]=0.0
           if(abs(kvec[n,2])<0.00000000000001):
             kvec[n,2]=0.0
           
  return kvec

def twist_input(filename,kpoints,index):
  with open(filename) as f:
    lines=f.readlines()

  with open("{}.{}".format(filename,index+1),'w') as f:
    i=0
    while(i < len(lines)):
      if lines[i].strip()[:8]=='K_POINTS':
        kgen_method=lines[i].split()[1]
        f.write("K_POINTS tpiba\n"); i+=1
        f.write("{}\n".format(kpoints.shape[0]))
        for n in range(kpoints.shape[0]):
          kvec=kpoints[n,:]
          f.write("{:06.16f} {:06.16f} {:06.16f} {:06.16f}\n".format(kvec[0],kvec[1],kvec[2],2.0/kpoints.shape[0]))
          if(kgen_method=='crystal'):
            i+=1
        if(kgen_method=='automatic'):
          i+=1
        break
      else:
        f.write(lines[i]); i+=1
          
  

def main():

    args=get_args()
    nk = np.array(args.kpoints)
    nt = np.array(args.twists)
    
    a,La = parse_input(args.file)

    Lb = reciprocal(La)

    Lb_sc = np.zeros((3,3))
    Lb_sc[0,:] = float(1.0/nk[0])*Lb[0,:]
    Lb_sc[1,:] = float(1.0/nk[1])*Lb[1,:]
    Lb_sc[2,:] = float(1.0/nk[2])*Lb[2,:]

    print()
    print('Lattice vectors (as rows):')
    print(La)
    print()
    print('Reciprocal lattice vectors: ')
    print(Lb)
    print()
    
    print("NK:")
    print(nk)

    kvec=generate_kpoint_grid(nk)
    print('The generated k-point grid:')
    kgrid=np.dot(kvec,Lb)
    print(kgrid)
    print()

    # tvec is to contain the actual twist grid, i.e. the grid
    # within the reciprocal unit cell of the supercell, while
    # tgrid will contain the twist grid as fractions of the lattice
    # vectors of the primitive cell. The latter is to be used for 
    # comparison with the SYMFILE vectors, which are in the same basis,
    # but produced for the twist generation.
    tvec=generate_kpoint_grid(nt)
    print('The generated twist grid (crystal coords.) :')
    for i in range(tvec.shape[0]):
      print("{}: {} {} {}".format(i,tvec[i,0],tvec[i,1],tvec[i,2]))
    print(tvec)

    tgrid=np.dot(tvec,Lb)
    tvec=np.dot(tvec,Lb_sc)
    
    print('The generated twist grid (cartesian coords.):')
    print(tvec)

    if(args.symfile=='nofile'):
      index_list=list(range(tgrid.shape[0]))
    else:
      print(args.symfile)
      sym=np.loadtxt(args.symfile)
      print(" ") 
      print("Unique twist vectors:")
      imax=sym.shape[0]
      ind=0
      index_list=[]
      for i in range(imax): # Loop over symmetry-reduced grid
        for j in range(tgrid.shape[0]): # Loop over full grid
          if(np.linalg.norm(sym[i,:3]-tgrid[j,:])<0.0000001):
            if(ind>imax-1):
              sys.exit("Found too many matches between full grid and symmetry-reduced grid.")
            print("tvec {}: ({:04.4f}, {:04.4f}, {:04.4f})  , weight= {}".format(j, tvec[j,0], tvec[j,1], tvec[j,2],sym[i,3]))
            ind+=1
            index_list.append(j)
      print("NOTE: This symmetry comparison works currently only for even k-grids, not for eg. 2x2x1")

    ind=0
    print("Generating the twisted input files.")
    for i in index_list:
      print(i)
      twist_input(args.file,kgrid+tvec[i,:],ind)
      ind+=1
    
    
    
if __name__ == '__main__':
    main()
