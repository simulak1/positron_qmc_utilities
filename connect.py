import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def get_args():
    """Define the task arguments with the default values.                         
    Returns:                                                                      
      experiment parameters                                                     
    """
    
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument(
        '--files',
        help='lineplot.dat-files.',
        type=str,
        nargs='+',
        default=['lineplot.dat']
    )
    
    # Data files arguments                                                        
    args_parser.add_argument(
        '--dft',
        help='Filename of the DFT result.',
        type=str,
        default='acar1d_100_ave'
    )

    args_parser.add_argument(
    '--weight',
    help=
        ''' 
        value by which the APMD coefficients are to be weighted.
        Hint: in si-units, S=1.9014321124319818e-21
        ''',
    type=float,
    default=1
  )

    args_parser.add_argument(
        '--ylog',
        help='log scale (1) or not (0).',
        type=int,
        default=0
    )

    args_parser.add_argument(
        '--ylim',
        help='ylim1,ylim2',
        type=str,
        default="not"
    )

    args_parser.add_argument(
        '--xlim',
        help='xlim1,xlim2',
        type=str,
        default="not"
    )

    
    return args_parser.parse_args()

def main():

    args=get_args()

    if(args.ylim!="not"):
        temp=args.ylim.split(',')
        if(len(temp)!=2):
            sys.exit("Wrong ylim input")
        ylim1=float(temp[0])
        ylim2=float(temp[1])

    if(args.xlim!="not"):
        temp=args.xlim.split(',')
        if(len(temp)!=2):
            sys.exit("Wrong ylim input")
        xlim1=float(temp[0])
        xlim2=float(temp[1])
        
    data=[]
    for fname in args.files:
        a=np.loadtxt(fname)
        print(a.shape)
        data.append(a[:,1])
    x=a[:,0]
    y=np.array(data)*args.weight
    print(y.shape)
    ymean=np.mean(y,axis=0)
    ystd=np.std(y,axis=0)
    
    plt.plot(x,ymean,'b-')
    plt.plot(x,ymean+ystd,'r--')
    plt.plot(x,ymean-ystd,'r--')
    
    plt.grid(True)
    plt.legend(["Mean data","60% confidence interval"])
    plt.xlabel("G-vector length (a.u.)")
    plt.ylabel("APMD value")
    if(args.ylim!="not"):
        plt.ylim([ylim1,ylim2])
    if(args.xlim!="not"):
        plt.xlim([xlim1,xlim2])
    if(not(args.ylog==0)):
        plt.yscale('log')
    plt.show()
    
if __name__ == '__main__':
    main()
