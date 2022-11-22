import numpy as np
import argparse
import matplotlib.pyplot as plt
import sys

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
    help='Number of twists in one direction (assumed all directions the same) for each file.',
    nargs='+',
    type=int,
    default=[2]
  )

  args_parser.add_argument(
    '--lineplots',
    help='Filenames of lineplot.dat.',
    nargs='+',
    type=str,
    default=[]
  )
  
  return args_parser.parse_args()


def parse(fn):
    with open(fn) as f:
        lines=f.readlines()
        for i in range(len(lines)):
            if(lines[i]=="Supercell reciprocal lattice vectors (au)\n"):
                pbmat=np.zeros((3,3))
                for j in range(3):
                    words=lines[i+j+1].split()
                    for k in range(3):
                        pbmat[j,k]=float(words[k])

            elif(lines[i]=="Number of G-vectors in set\n"):
                Ng=int(lines[i+1])
            elif(lines[i]=="G-vector components Gx, Gy, Gz (au)\n"):
                gvec=np.zeros((Ng,3))
                for j in range(Ng):
                    words=lines[i+j+1].split()
                    for k in range(3):
                        gvec[j,k]=float(words[k])
            elif(lines[i]=="APMD coefficients (mean value, standard error)\n"):
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
  g=[]; p=[]; e=[]
  g.append(np.linalg.norm(gvec[0,:]))
  p.append(pmd[0,0])
  e.append(pmd[0,1])
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
      e.append(pmd[i,1])

  g=np.array(g)
  p=np.array(p)
  e=np.array(e)
  isort=np.argsort(g)

  return g[isort],p[isort],e[isort]

def read_lineplot(file):
  x=[]; y=[]
  with open(file) as f:
    lines=f.readlines()
    for line in lines:
      words=line.split()
      x.append(float(words[0]))
      y.append(float(words[1]))

  return x,y
      
def main():

    args=get_args()
    temp=args.angle.split(',')
    ia=int(temp[0]); ib=int(temp[1]); ic=int(temp[2])

    nt=args.ntwist
    if(len(nt)!=len(args.files)):
       sys.exit("Number of twists must be provided for all files")
    
    dict={}

    ind=0
    plt.figure(1)
    mlist=['r+','bo','g*','c.','m<','y>','wh','kx','bd']
    for fname in args.files:
      pbmat,Ng,gvec,pmd=parse(fname)
      #pmd=pmd/np.max(pmd)

      A=np.transpose(pbmat)/float(nt[ind])
      AtA=np.matmul(np.transpose(A),A)
      mpinv=np.matmul(np.linalg.inv(AtA),np.transpose(A))

      g,p,e=linepoints(gvec,pmd,mpinv,ia,ib,ic)
      p=p/np.max(p)
      print("G-vector, value, error:")
      for i in range(g.shape[0]):
        print("{:04.3f}, {:.3e}, {:.3e}".format(g[i],p[i],e[i]))
      plt.errorbar(g,p,e,color='green',marker='+')#,mlist[ind])
      ind+=1

    lines=['r-','b-','g-']
    ind=0
    for fname in args.lineplots:
      x,y=read_lineplot(fname)
      y=y/np.max(y)
      x=np.array(x); y=np.array(y); #y=y/np.max(y)
      plt.plot(x,y,lines[ind])
      ind+=1
      
    plt.title("Points found in direction "+args.angle)
    plt.grid()
    plt.legend(args.files)
    plt.show()

if __name__=='__main__':
    main()
