import numpy as np
import argparse
import matplotlib.pyplot as plt
import sys

'''
This script searches for g-vectors along a given line, takes APMD values
corresponding to the vectors from one or multiple expval.data-files, and 
plots the values as a function of the g-vector lentghts. If multiple 
expval.data-files are fed in to the script, the plto will also contain 
errorbars of the APMD values. 
'''

def get_args():

  args_parser = argparse.ArgumentParser()

  args_parser.add_argument(
    '--files',
    help='Filenames.',
    nargs='+',
    type=str,
    default=['expval.data']
  )

  args_parser.add_argument(
    '--angle',
    help='3 integers: i,j,k, so that we look points in the direction i*b1+j*b2+k*b3.',
    type=str,
    default='1,0,0'
  )

  args_parser.add_argument(
    '--ntwist',
    help='Number of twists at each lattice vector direction, given as i,j,k.',
    type=str,
    default='1,1,1'
  )

  args_parser.add_argument(
    '--weight',
    help='''
    value by which the APMD coefficients are to be weighted.
    Hint: in si-units, S=1.9014321124319818e-21
''',
    type=float,
    default=1
  )
  
  return args_parser.parse_args()


def parse(fn,nt):
    '''
    Parses expval.data-file and returns the reciprocal lattice vectors,
    Number of g-vectors, list of g-vectors and the APMD values. Because
    we also want to parse merged expval.data-files, where multiple twisted
    g-vector-grids are merged, we give as an argument to this function the 
    array nt, that describes the number of twists at each lattice vector
    direction. Then the i:th reciprocal lattice vector is divided by nt[i],
    so that by integer linear combinations of the newly generated reciprocal
    lattice vectors we can found ALL of the g-vectors, including the twisted 
    ones. 
    '''
    with open(fn) as f:
        lines=f.readlines()
        for i in range(len(lines)):
            if(lines[i]=="Supercell reciprocal lattice vectors (au)\n"):
                pbmat=np.zeros((3,3))
                for j in range(3):
                    words=lines[i+j+1].split()
                    for k in range(3):
                        pbmat[j,k]=float(words[k])/nt[j]

            elif(lines[i]=="Number of G-vectors in set\n"):
                Ng=int(lines[i+1])
            elif(lines[i]=="G-vector components Gx, Gy, Gz (au)\n"):
                gvec=np.zeros((Ng,3))
                for j in range(Ng):
                    words=lines[i+j+1].split()
                    for k in range(3):
                        gvec[j,k]=float(words[k])
            elif(lines[i]=="Complex pair-density coefficients (real part, imaginary part)\n"):
                pmd=np.zeros((Ng,2))
                for j in range(Ng):
                    words=lines[i+j+1].split()
                    for	k in range(2):
                        pmd[j,k]=float(words[k])

    return pbmat,Ng,gvec,pmd

def testint(r):
  tol=0.0000001
  if(abs(round(r[0])-r[0])>tol):
    return -1
  elif(abs(round(r[1])-r[1])>tol):
    return -1
  elif(abs(round(r[2])-r[2])>tol):
    return -1
  else:
    return 1

def linepoints(gvec,pmd,A,ia,ib,ic):
  g=[]; p=[]
  g.append(np.linalg.norm(gvec[0,:]))
  p.append(pmd[0,0])
  for i in range(1,gvec.shape[0]):
    gv=gvec[i,:]
    intvec=np.matmul(A,gv)
    if(testint(intvec)<0):
      sys.exit("Non-integer linear combination of reciprocal lattice vectors:[{} {} {}]".format(intvec[0],intvec[1],intvec[2]))
    ivec=np.array((round(intvec[0]),round(intvec[1]),round(intvec[2])))
    absivec=np.absolute(ivec)
    if(np.nonzero(absivec)[0].shape[0]>1):
      ii=np.min(absivec[np.nonzero(absivec)])
    else:
      ii=absivec[np.nonzero(absivec)]

    intvec=intvec/ii
    if(testint(intvec)<0):
      continue
    if((round(intvec[0])==ia and round(intvec[1])==ib and round(intvec[2])==ic)):
      g.append(np.sqrt(gv[0]**2+gv[1]**2+gv[2]**2))
      p.append(pmd[i,0])

  g=np.array(g)
  p=np.array(p)
  
  isort=np.argsort(g)

  return g[isort],p[isort]

def main():

    args=get_args()
    print(" ")
    print("====== Plot APMD values along a given line =====")
    print(" ")
    print("CASINO expval files to be parsed:")
    for file in args.files:
      print(file)
    print(" ")

    temp=args.angle.split(',')
    ia=int(temp[0]); ib=int(temp[1]); ic=int(temp[2])

    temp=args.ntwist.split(',')
    nt1=int(temp[0]); nt2=int(temp[1]); nt3=int(temp[2])
    nt=np.array([nt1,nt2,nt3])
    
    dict={}

    ind=0
    plt.figure(1)
    pointlists=[]
    for fname in args.files:
      print("Parsing {}...".format(fname))
      pbmat,Ng,gvec,pmd=parse(fname,nt)
      
      A=np.transpose(pbmat)
      AtA=np.matmul(np.transpose(A),A)
      mpinv=np.matmul(np.linalg.inv(AtA),np.transpose(A))

      g,p=linepoints(gvec,pmd,mpinv,ia,ib,ic)
      p=p*args.weight
      
      pointlists.append(p)
      ind+=1
      

    for p in pointlists[1:]:
      if p.shape[0]!=pointlists[0].shape[0]:
        sys.exit("Error: Different number of points along a line found from expval.data-files. ")
      
    YX=np.array(pointlists)

    yval=np.mean(YX,axis=0)
    yvar=np.std(YX,axis=0)
    print(" ")
    print("Obtained statistics: ")
    print("g-vec length, mean value, error")
    for i in range(yval.shape[0]):
          print("{}, {}, {}".format(g[i],yval[i],yvar[i]))
    plt.errorbar(g,yval,yvar)
    plt.title("Points found in direction "+args.angle)
    plt.grid()
    plt.legend(args.files)
    plt.show()

if __name__=='__main__':
    main()
