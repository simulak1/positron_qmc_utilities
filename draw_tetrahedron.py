# -*- coding: utf-8 -*-
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import argparse

def get_args():
  """Define the task arguments with the default values.                         
  Returns:                                                                      
      experiment parameters                                                     
  """

  args_parser = argparse.ArgumentParser()

  # Data files arguments                                                        
  args_parser.add_argument(
    '--vfile',
    help='Local paths to file of vertices.',
    default='vertices.dat'
  )

  args_parser.add_argument(
    '--tfile',
    help='Local paths to file of tetrahedra.',
    default='tetrahedra.dat'
  )

  args_parser.add_argument(
    '--point',
    help='Plot tetrahedra related to a point',
    type=int,
    default=0
  )  

  return args_parser.parse_args()


def _read_vertices(filename='fort.66'):
    vertices=[]
    with open(filename) as f:
        lines=f.readlines()
        for line in lines:
            words=line.split()
            vertex=[float(words[0]),float(words[1]),float(words[2])]
            vertices.append(np.array(vertex))
    return np.array(vertices)

def _read_tetrahedra(filename='fort.76'):
    tetrahedra=[]
    with open(filename) as f:
        lines=f.readlines()
        for line in lines:
            words=line.split()
            tetrahedron=[int(words[0]),int(words[1]),int(words[2]),int(words[3])]
            tetrahedra.append(np.array(tetrahedron))
    return np.array(tetrahedra)

def plot_vertices(v,vind,ax):

    ax.scatter(v[vind-1,0],v[vind-1,1],v[vind-1,2],'ko',s=40)
    for i in vind:
        
        ax.text(v[i-1,0],v[i-1,1],v[i-1,2],  '%s' % (str(i)))

def plot_tetrahedron(t,v,ax):

    color=[np.random.rand(),np.random.rand(),np.random.rand(),1]
    #ab    
    ax.plot([v[t[0],0],v[t[1],0]],[v[t[0],1],v[t[1],1]],[v[t[0],2],v[t[1],2]],'-',color=color)
    #ac    
    ax.plot([v[t[0],0],v[t[2],0]],[v[t[0],1],v[t[2],1]],[v[t[0],2],v[t[2],2]],'-',color=color)
    #ad    
    ax.plot([v[t[0],0],v[t[3],0]],[v[t[0],1],v[t[3],1]],[v[t[0],2],v[t[3],2]],'-',color=color)
    #bc    
    ax.plot([v[t[1],0],v[t[2],0]],[v[t[1],1],v[t[2],1]],[v[t[1],2],v[t[2],2]],'-',color=color)
    #cd    
    ax.plot([v[t[2],0],v[t[3],0]],[v[t[2],1],v[t[3],1]],[v[t[2],2],v[t[3],2]],'-',color=color)
    #db    
    ax.plot([v[t[3],0],v[t[1],0]],[v[t[3],1],v[t[1],1]],[v[t[3],2],v[t[1],2]],'-',color=color)    
    
def get_points_in_complex(v,t):
    points=[]
    for i in range(t.shape[0]):
        for j in range(4):
            point_candidate = t[i,j]
            if(not(point_candidate in points) and not(point_candidate==-1)):
                points.append(point_candidate)
                
    return np.array(points)
    
def main():
    
    args=get_args()
    
    point=args.point
        
    v=_read_vertices(args.vfile)  
    t=_read_tetrahedra(args.tfile)    
    vind=get_points_in_complex(v,t)
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if(point==0):
      plot_vertices(v,vind,ax)
        
    print('file: {}'.format(args.vfile))
    for i in range(0,t.shape[0]):
      if(point>0):
        if(t[i,-1]!=-1 and (t[i,0]==point or t[i,1]==point or t[i,2]==point or t[i,3]==point)):
          print('Tetrahedron {}: {},{},{},{}'.format(i,t[i,0],t[i,1],t[i,2],t[i,3]))
          plot_tetrahedron(t[i,:]-1,v,ax)
          plot_vertices(v,t[i,:],ax)
      else:
        print('Tetrahedron {}: {},{},{},{}'.format(i,t[i,0],t[i,1],t[i,2],t[i,3]))
        if(t[i,-1]!=-1):
          plot_tetrahedron(t[i,:]-1,v,ax)
                
    plt.show()

if __name__ == '__main__':
    main()


