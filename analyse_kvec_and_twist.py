from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import re

def get_args():
  """Define the task arguments with the default values.                         
  Returns:                                                                      
      experiment parameters                                                     
  """

  args_parser = argparse.ArgumentParser()

  # Data files arguments                                                        
  args_parser.add_argument(
      '--files',
      help='Filenames.',
      nargs='+',
      type=str,
      default=['k_offsets']
  )

  args_parser.add_argument(
    '--invert2direct',
    help='''
    In case the k_offsets-files are in cartesian coordinates, invert them
    to a coordinate system with a reciprocal lattice vector basis. Requires
    The reciprocal lattice vectors.
    ''',
    type=int,
    default=0
  )

  args_parser.add_argument(
    '--latvec',
    help='Lattice vector components in the following syntax_: a1,a2,a3,b1,b2,b3,c1,c2,c3.',
    type=str,
    default='1,0,0,0,1,0,0,0,1'
  )

  args_parser.add_argument(
    '--label',
    help='give indices as a comma-separated list of points that you would like to label in plot.',
    type=str,
    default='-1'
  )

  args_parser.add_argument(
    '--labelfile',
    help='Name of the file that needs to be labeled.',
    type=str,
    default='k_offset'
  )
  
  return args_parser.parse_args()

def parse(file):
    
    with open(file) as f:
        lines=f.readlines()
        kvec=np.zeros((len(lines)-1,3))
        i=0
        Nk=int(lines[0])
        for line in lines[1:]:
            words=line.split()
            kvec[i,:]=np.array((float(words[0]),float(words[1]),float(words[2])))
            i+=1
            if(i>Nk-1):
              break
    return kvec

def invert2direct(args,grids):
  latvec=args.latvec.split(',')
  if len(latvec)!=9:
    sys.exit('Wrong number of lattice vector components specified. ')
  amat=np.zeros((3,3))
  for i in range(9):
    ix=i%3
    iy=int((i-i%3)/3)
    amat[iy,ix]=float(latvec[i]) # NOTE!                                                              

  AtA=np.matmul(np.transpose(amat),amat)
  invAtA=np.linalg.inv(AtA)
  for file in args.files:
    grids[file]=np.transpose(np.matmul(np.matmul(invAtA,amat),np.transpose(grids[file])))
  return grids
  
def main():
  
  args=get_args()
    
  grids={}

  for file in args.files:
    grids[file]=parse(file)
      
  if(args.invert2direct==1):
    grids=invert2direct(args,grids)

  tl=args.label.split(',')
  labels=[]
  for i in range(len(tl)):
    labels.append(int(tl[i])-1)
    
  # ---------------- PLOTS --------------------#
  fig1 = plt.figure(4)
  ax1 = fig1.add_subplot(111, projection='3d')
  for file in args.files:
    ax1.scatter(grids[file][:,0],grids[file][:,1],grids[file][:,2],alpha=0.5)
    if labels[0]>0  and args.labelfile==file:
      ax1.scatter(grids[file][labels,0],grids[file][labels,1],grids[file][labels,2],'k*',linewidth=5)
      for il in labels:
        ax1.text(grids[file][il,0],grids[file][il,1],grids[file][il,2],  '%s' % (str(il+1)))
    ax1.legend(args.files)
  #plt.show()

  fig,ax=plt.subplots(3,1)
  for i in range(3):
    if(i==0):
      ix=0;iy=1;title='xy-plot'
    elif(i==1):
      ix=0;iy=2;title='xz-plot'
    elif(i==1):
      ix=1;iy=2;title='yz-plot'
    for file in args.files:
      ax[i].scatter(grids[file][:,ix],grids[file][:,iy])
      #ax[i].legend(args.files)
      ax[i].set_title(title)
      ax[i].grid(True)
  plt.show()

    
    

if __name__=='__main__':
    main()
