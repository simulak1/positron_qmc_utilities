import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def get_args():
    """
    Define the task arguments with the default values.
    Returns: experiment parameters     
    """
    
    args_parser = argparse.ArgumentParser()
    
    # Data files arguments

    args_parser.add_argument(
        '--fdoppler',
        help='Filename of the doppler measurement.',
        type=str,
        default='none'
    )

    args_parser.add_argument(
        '--dft',
        help='Filename of the DFT result.',
        nargs='+',
        type=str,
        default=[]
    )

    args_parser.add_argument(
        '--qmc',
        help='Filename of the qmc result.',
        nargs='+',
	type=str,
        default=['lineplot.dat']
    )

    args_parser.add_argument(
        '--fcore',
        help='Filename of the core electron APMD.',
        type=str,
        nargs='+',
        default=['no']
    )

    args_parser.add_argument(
	'--plotcore',
        help='Plot core spectrum separately.',
	type=int,
	default=-1
    )

    args_parser.add_argument(
        '--plotsw',
        help='Plot limits of the S- and W-parameters.',
	type=int,
        default=-1
    )
    
    args_parser.add_argument(
        '--chwidth',
        help='Channel width  of the experimental data (in au.)',
        type=float,
        default=1
    )
    
    args_parser.add_argument(
        '--vnorm',
        help='The valence electron contribution.',
        type=float,
        default=1
    )

    args_parser.add_argument(
        '--fwhm',
        help='Full width at half maximum (keV) of the resolution function, for convolution.',
	type=float,
        default=0.85
    )

    args_parser.add_argument(
        '--convrange',
        help='The width of the convolution window: [-convrange,convrange].',
	type=float,
        default=5.0
    )

    args_parser.add_argument(
        '--maxmom',
	help='The maximum momentum value included in the projections.',
        type=float,
        default=4.0
    )
    
    args_parser.add_argument(
        '--ylog',
        help='Plot in log scale (1) or not (0).',
        type=int,
        default=0
    )


    args_parser.add_argument(
        '--errorplot',
        help='are errors contained in the lineplot.dat-file? 1=yes, 0=no.',
	type=int,
        default=1
    )

    args_parser.add_argument(
        '--angle',
	help='Name of the projection angle.',
        type=str,
        default='?'
    )

    args_parser.add_argument(
        '--srange',
        help='The window for estimating S-parameter, given as a string a,b in eV.',
	type=str,
        default='0.0,0.4'
    )

    args_parser.add_argument(
        '--wrange',
        help='The window for estimating W-parameter, given as a string a,b in eV.',
	type=str,
        default='1.6,2.4'
    )

    args_parser.add_argument(
        '--figname',
        help='If given, the image to be plotted will be saved under this name.',
        type=str,
        default=''
    )

    args_parser.add_argument(
        '--skipconv',
        help='obvious: >1=yes, <0=no, default=no',
        type=int,
        default=-1
    )
    
    return args_parser.parse_args()

def get_imax(r,rmax):
    imax=r.shape[0]
    for i in range(len(r)):
        if(r[i]>rmax):
            imax=i
            return imax
    return imax

def convolute(x,y,fwhm,convrange,skip=-1):
    if(skip<0):
        # x and y are in range [0,xmax]. Mirror it now to range [-xmax,xmax].
        y2=np.concatenate((np.flip(y),y))
        x2=np.concatenate((-np.flip(x),x))

        # Compute the convolution function
        xx=np.arange(-convrange,convrange,x[1]-x[0])
        if(convrange>np.max(x)):
            sys.exit("Pick a smaller convolution range.")
        sigma=fwhm/2.355
        convf=np.exp(-1.0*(xx**2)/(2*sigma**2))
    
        # Convolve
        convoluted=np.convolve(y2,convf,mode='same')
        
        # Remove the mirror image
        y2=convoluted[-x.shape[0]:]
    
        # Normalize
        area=np.trapz(y2,x)
        y2=y2/area

        return np.array([x,y2])
    else:
        return np.array([x,y])

def load_data(args,filelist,is_core,vnorm,core_data,sr=0):
    if(is_core and len(core_data[0])==1 and core_data[0][0]==0):
        sys.exit("Error: Trying to add core data to projection, but there is none.")
    lp=[]
    index=0
    for file in filelist:
        data=np.transpose(np.loadtxt(file,skiprows=sr))
        data=data[:,:get_imax(data[0,:],args.maxmom)]
        norm=vnorm/np.trapz(data[1,:],x=data[0,:])
        data[1,:]*=norm
        if(data.shape[0]==3):
            data[2,:]*=norm
        if(is_core):
            core=np.interp(data[0,:],core_data[index][:,0],core_data[index][:,1])
            data[1,:]+=core
            index+=1
        lp.append(data)
        
    data=np.array(lp)
    return data

def plot_data(args,i,dft,qmc1,qmc2,qmc3,core_data,sa,sb,wa,wb,ydop=0,xx=0):
    # qmc1= mean value, qmc2=minimum, qmc3=maximum
    ax1=plt.subplot(212)
    #ax1=plt.subplot(111)
    ax2 = plt.subplot(221)
    ax3 = plt.subplot(222)
    legends=[]

    # Set ylimits for S- and W-parameter windows
    sind1=get_imax(qmc1[0,:],sa); sind2=get_imax(qmc1[0,:],sb);
    sylim1=1.07*qmc1[1,sind1:sind2].max()
    sylim2=.93*qmc1[1,sind1:sind2].min()

    wind1=get_imax(qmc1[0,:],wa); wind2=get_imax(qmc1[0,:],wb);
    wylim1=1.1*qmc1[1,wind1:wind2].max()
    wylim2=.9*qmc1[1,wind1:wind2].min()
    
    for j in range(len(dft)):
        ax1.plot(dft[j][0,:],dft[j][1,:],'--')
        ax2.plot(dft[j][0,:],dft[j][1,:],'--')
        ax3.plot(dft[j][0,:],dft[j][1,:],'--')
        legends.append(args.dft[j].split("/")[-1])
    ax1.set_xlim([0,2.5])
    ax2.set_xlim([sa,sb])
    ax2.set_ylim([sylim2,sylim1])
    ax3.set_xlim([wa,wb])
    ax3.set_ylim([wylim2,wylim1])
    ax1.plot(qmc1[0,:],qmc1[1,:],'-')
    ax2.plot(qmc1[0,:],qmc1[1,:],'-')
    ax3.plot(qmc1[0,:],qmc1[1,:],'-')
    ax1.fill_between(qmc1[0,:],qmc2[1,:],qmc3[1,:],color='gray', alpha=0.2)
    ax2.fill_between(qmc1[0,:],qmc2[1,:],qmc3[1,:],color='gray', alpha=0.2)
    ax3.fill_between(qmc1[0,:],qmc2[1,:],qmc3[1,:],color='gray', alpha=0.2)

    legends.append(args.qmc[i].split("/")[-1])
    
    if(args.plotcore>0):
        for i in range(len(core_data)):
            imax=get_imax(core_data[i][:,0],args.maxmom)
            ax1.plot(core_data[i][:imax,0],core_data[i][:imax,1],'r')
            ax3.plot(core_data[i][:,0],core_data[i][:,1],'r')
            legends.append('Core spectrum')

    if(args.fdoppler!='none'):
        ax1.plot(xx,ydop,'k-+')
        ax2.plot(xx,ydop,'k-+')
        ax3.plot(xx,ydop,'k-+')
        legends.append('Experiment')

    if(args.plotsw>0):
        ax1.axvline(sb,ls='--')
        ax1.axvline(wa,ls='--')
        ax1.axvline(wb,ls='--')
        
    if(not(args.ylog==0)):
        ax1.set_yscale('log')
    ax1.grid()
    ax1.set_ylim([1.*10**-4,sylim1])
    ax1.set_title('Diamond pmd data projection to {}'.format(args.angle))
    ax2.set_title("S-parameter area")
    ax3.set_title("W-parameter area")
    
    ax1.set_xlabel('Momentum (a.u.)')
    ax1.legend(legends)
    if(len(args.figname)>0):
        plt.savefig(args.figname, dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', format='pdf',
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    metadata=None)

def plot_data2(args,i,dft,qmc1,qmc2,qmc3,core_data,sa,sb,wa,wb,ydop=0,xx=0):

    def _discrepancy(a,b,x):
        newa=[]
        for i in range(x.shape[0]):
            ind=np.argmin(np.abs(a[0,:]-x[i]))
            newa.append(abs(a[1,ind]-ydop[i]))
        return np.array(newa)
            
    ax1=plt.subplot(111)
    legends=[]

    sind1=get_imax(qmc1[0,:],sa); sind2=get_imax(qmc1[0,:],sb);
    sylim1=1.07*qmc1[1,sind1:sind2].max()
    for j in range(len(dft)):
        ddft=_discrepancy(dft[j],ydop,xx)
        ax1.plot(xx,ddft,'-')
        legends.append(args.dft[j].split("/")[-1])
    ax1.set_xlim([0,2.5])
    dqmc=_discrepancy(qmc1,ydop,xx)
    ax1.plot(xx,dqmc,'-')
    #ax1.fill_between(qmc1[0,:],qmc2[1,:],qmc3[1,:],color='gray', alpha=0.2)

    
    legends.append(args.qmc[i].split("/")[-1])
    
    if(not(args.ylog==0)):
        ax1.set_yscale('log')
    ax1.grid()
    #ax1.set_ylim([1.*10**-4,sylim1])
    ax1.set_title('Diamond pmd data projection to {}'.format(args.angle))
    
    ax1.set_xlabel('Momentum (a.u.)')
    ax1.legend(legends)
    if(len(args.figname)>0):
        plt.savefig(args.figname, dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', format='pdf',
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    metadata=None)


def plot_data3(args,i,dft,qmc1,qmc2,qmc3,core_data,sa,sb,wa,wb,ydop=0,xx=0):

    def _discrepancy(a,b,x):
        newa=[]
        for i in range(x.shape[0]):
            ind=np.argmin(np.abs(a[0,:]-x[i]))
            newa.append(a[1,ind]-ydop[i])
        newa=np.array(newa)
        mse=np.sum((newa-np.mean(newa))**2)/newa.shape[0]
        return newa,mse
    
    f, (ax1, ax2) = plt.subplots(2, 1,
                                 gridspec_kw={'height_ratios': [3, 1]},
                                 sharex=True)
    #f.tight_layout()
    plt.subplots_adjust(hspace=.05,wspace=1)
    legends=[]; legends2=[]

    
    # Set ylimits for S- and W-parameter windows
    sind1=get_imax(qmc1[0,:],sa); sind2=get_imax(qmc1[0,:],sb);
    sylim1=1.07*qmc1[1,sind1:sind2].max()
    sylim2=.93*qmc1[1,sind1:sind2].min()

    wind1=get_imax(qmc1[0,:],wa); wind2=get_imax(qmc1[0,:],wb);
    wylim1=1.1*qmc1[1,wind1:wind2].max()
    wylim2=.9*qmc1[1,wind1:wind2].min()
    dft_names=["State-dependent DFT","State-independent DFT"]
    for j in range(len(dft)):
        ddft,mse_dft=_discrepancy(dft[j],ydop,xx)
        ax2.plot(xx,ddft,'-')
        ax1.plot(dft[j][0,:],dft[j][1,:],'--')
        #legends.append(args.dft[j].split("/")[-1])
        legends.append(dft_names[j])
        legends2.append("MSE: {:.4e}".format(mse_dft))
    ax1.set_xlim([0,2.5])
    dqmc,mse_qmc=_discrepancy(qmc1,ydop,xx)
    dqmc2,mse_qmc2=_discrepancy(qmc2,ydop,xx)
    dqmc3,mse_qmc3=_discrepancy(qmc3,ydop,xx)
    legends2.append("MSE: {:.4e}".format(mse_qmc))
    ax2.plot(xx,dqmc,'-')
    ax1.plot(qmc1[0,:],qmc1[1,:],'-')
    ax1.fill_between(qmc1[0,:],qmc2[1,:],qmc3[1,:],color='gray', alpha=0.2)
    ax2.fill_between(xx,dqmc2,dqmc3,color='gray', alpha=0.2)

    legends.append("QMC")
    
    for i in range(len(core_data)):
            imax=get_imax(core_data[i][:,0],args.maxmom)
            ax1.plot(core_data[i][:imax,0],core_data[i][:imax,1],'r')
            #ax2.plot(core_data[i][:,0],core_data[i][:,1],'r')
            legends.append('Core spectrum')

    if(args.fdoppler!='none'):
        ax1.plot(xx,ydop,'k-+')
        legends.append('Experiment')

    if(args.plotsw>0):
        ax1.axvline(sb,ls='--')
        ax1.axvline(wa,ls='--')
        ax1.axvline(wb,ls='--')
        
    if(not(args.ylog==0)):
        ax1.set_yscale('log')
        ax2.set_yscale('log')

    ax1.grid()
    ax1.set_ylabel("APMD projection")
    ax2.grid()
    ax2.set_ylabel("Error to experiment")
    ax1.set_ylim([1.*10**-4,sylim1])
    #ax1.set_title('Momentum projection')
    
    
    ax2.set_xlabel('Momentum (a.u.)')
    ax1.legend(legends,)
    ax2.legend(legends2,loc="upper right")
    if(len(args.figname)>0):
        plt.savefig(args.figname, dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', format='pdf',
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    metadata=None)


        
def getSW(X,a,b):
    from scipy.integrate import simps

    ai=np.argmin(np.absolute(X[0,:]-a))
    bi=np.argmin(np.absolute(X[0,:]-b))

    return simps(X[1,ai:bi+1],X[0,ai:bi+1])
    
def main():

    args=get_args()
    
    # fwhm for convolution, keV transformed to... a.u.
    fwhm=0.54*args.fwhm
    # Range of the convolution, keV, transformed...
    grange=0.54*args.convrange

    # S-parameter window
    tempstr=args.srange.split(',')
    if(len(tempstr)<2 or len(tempstr)>2):
        sys.exit("Error: check the format of the srange-parameter given in command line.")
    sa=float(tempstr[0])
    sb=float(tempstr[1])

    # W-parameter window
    tempstr=args.wrange.split(',')
    if(len(tempstr)<2 or len(tempstr)>2):
        sys.exit("Error: check the format of the wrange-parameter given in command line.")
    wa=float(tempstr[0])
    wb=float(tempstr[1])
    
    plot_error=False
    if(args.errorplot==1):
        plot_error=True
    
    is_core=not(args.fcore[0]=='no')
    core_data=[]
    if(is_core):
        for name in args.fcore:
            core_data.append(np.loadtxt(name,skiprows=1))
    else:
        core_data.append(np.array((0,)))
        
    if(args.plotcore>0):
        if(not(is_core)):
            sys.exit("Required core spectrum plot, but no core file provided.")
        
    dft_data=load_data(args,args.dft,is_core,args.vnorm,core_data=core_data,sr=2)
    qmc_data=load_data(args,args.qmc,is_core,args.vnorm,core_data=core_data)
    # Load, cut and normalize the experimental data, if provided.
    have_experiment=False
    if(args.fdoppler!='none'):
        have_experiment=True
        ydop=np.loadtxt(args.fdoppler)
        if(len(ydop.shape)>1):
            print("WARNING: Channels are included in doppler data. Check how you treat them!")
            xx=ydop[:,0]/3.91*0.54
            ydop=ydop[:,1]
        else:
            xx=np.arange(0,7,args.chwidth)
        xx=xx[:get_imax(xx,args.maxmom)]
        ydop2=ydop[:xx.shape[0]]
        norm=1./np.trapz(ydop2,x=xx)
        ydop2*=norm

        print(" ")
        print("*** Experimental S- and W-parameters ***")
        print(" ")
        print("- file: {}".format(args.fdoppler))
        print("S: {:04.4f}".format(getSW(np.array([xx,ydop2]),sa,sb)))
        print("W: {:04.4f}".format(getSW(np.array([xx,ydop2]),wa,wb)))
    
        
    print(" ")
    print("*** S- and W-parameters (DFT) ***")
    print(" ")
    dftc=[]
    
    for i in range(len(args.dft)):
        print("- file: {}".format(args.dft[i]))
        dftconv=convolute(dft_data[i,0,:],dft_data[i,1,:],fwhm,grange,args.skipconv)
        dftc.append(dftconv)
        print("S: {:04.4f}".format(getSW(dftconv,sa,sb)))
        print("W: {:04.4f}".format(getSW(dftconv,wa,wb)))

    print(" ")
    print("*** S- and W-parameters (QMC) ***")
    print(" ")
    for i in range(len(args.qmc)):
        qmcconv1=convolute(qmc_data[i,0,:],qmc_data[i,1,:],fwhm,grange,args.skipconv)
        qmcconv2=convolute(qmc_data[i,0,:],qmc_data[i,1,:]-qmc_data[i,2,:],fwhm,grange,args.skipconv)
        qmcconv3=convolute(qmc_data[i,0,:],qmc_data[i,1,:]+qmc_data[i,2,:],fwhm,grange,args.skipconv)
        S1=getSW(qmcconv1,sa,sb)
        S2=getSW(qmcconv2,sa,sb)
        S3=getSW(qmcconv3,sa,sb)
        W1=getSW(qmcconv1,wa,wb)
        W2=getSW(qmcconv2,wa,wb)
        W3=getSW(qmcconv3,wa,wb)
        print("Angle {}: S= {:04.4f} +- {:04.4f}".format(args.angle,S1,(S2-S3)/2.0))
        print("Angle {}: W= {:04.4f} +- {:04.4f}".format(args.angle,W1,(W2-W3)/2.0))
        if(have_experiment):
            plot_data(args,i,dftc,qmcconv1,qmcconv2,qmcconv3,core_data,sa,sb,wa,wb,ydop2,xx)
        else:
            plot_data(args,i,dftc,qmcconv1,qmcconv2,qmcconv3,core_data,sa,sb,wa,wb)
    print(" ")
    if(len(args.figname)<1):
        plt.show()
    
if __name__ == '__main__':
    main()
