"""
python welltdf.py BO-02_Sonic.ascii --dtcols 0 1 --slowness --headerlines 28
python welltdf.py BO-02_Sonic.ascii --dtcols 0 1 --slowness --headerlines 28 --zrange 0 1200 2200
python welltdf.py BO-02_TD.ascii --dtcols 2 3  --headerlines 14 --dtmultiplier -1 -0.0005
python welltdf.py BO-02_TD.ascii --dtcols 2 3  --headerlines 14 --dtmultiplier -1 -0.0005 --hideplot
python welltdf.py BO-02_TD.ascii --dtcols 2 3 --headerlines 14 --dtmultiplier -1 -0.0005 --zrange 1 0 1510
python welltdf.py BO-02_TD.ascii --dtcols 2 3 --headerlines 14 --dtmultiplier -1 -0.0005 --qcfcorrection 1 1.32766637  0.9283177

"""

import os.path
import datetime
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import scipy.stats as st
import pandas as pd

def topsin(fname,zc,topc):
    """."""
    topz = np.genfromtxt(fname,usecols=zc)
    tops = np.genfromtxt(fname,usecols=topc,dtype=None)
    return topz,tops

def listtops(topz,tops):
    """."""
    for i in range(topz.size):
        print('%10.0f  %20s' % (topz[i],tops[i]))
    return

def listitops(topz,topcsum,wname,tops,xc,yc):
    """."""
    for i in range(topz.size):
        vav = topz[i] / topcsum[i]
        print('%12.2f  %12.2f  %10.0f  %10.3f   %10.0f   %10.0f   %20s  %20s' %\
        (xc,yc,topz[i],topcsum[i],topcsum[i] *2000.0,vav,wname,tops[i]))
    return

def zprsin(fname,nh,zc,tpc,rhoc,tsc):
    """."""
    zprs = np.genfromtxt(fname,usecols=(zc,tpc,rhoc,tsc),skip_header=nh)
    return zprs[:,0],zprs[:,1],zprs[:,2],zprs[:,3]

def zprin(fname,nh,zc,tpc,rhoc):
    """."""
    zpr = np.genfromtxt(fname,usecols=(zc,tpc,rhoc),skip_header=nh)
    return zpr[:,0],zpr[:,1],zpr[:,2]

def zpsin(fname,nh,zc,tpc,tsc):
    """."""
    zps = np.genfromtxt(fname,usecols=(zc,tpc,tsc),skip_header=nh)
    return zps[:,0],zps[:,1],zps[:,2]


def zpin(fname,nh,zc,tpc,nulval):
    """."""
    zp = np.genfromtxt(fname,usecols=(zc,tpc),skip_header=nh)
    return zp[:,0],zp[:,1]

def listdata(z,t,vi):
    """."""
    for i in range(len(z)):
        print('%10.3f  %10.6f  %10.0f' %(z[i],t[i],vi[i]))
    return

def listindata(z,t):
    """."""
    for i in range(len(z)):
        print('%10.1f  %10.3f ' % (z[i],t[i]))
    return

def listheader(fname):
    """."""
    lasf = open(fname,'r')
    nh = 0
    for line in lasf:
        if line[0] == '~' and line[1] == 'A':
            print(line[:-1])
            nh += 1
            return nh
        else:
            print(line[:-1])
            nh += 1
    return

def list_dt(z,t):
    """."""
    for i in range(z.size):
        print('%10.2f  %10.3f' % (z[i],t[i]))
    return

def lwzfit(z,vi,topsz,topsname):
    """."""
    f0 = open('lwzcoef.lst','w')
    now = datetime.datetime.now()
    print("#", now.strftime("%Y-%b-%d  %H:%M"), file=f0)
    print("# Intercept  Slope  r_value  p_value  std_err topname min_depth  max_depth", file=f0)
    slope = np.zeros(z.size)
    intercept = np.zeros(topsz.size - 1)
    r_value = np.zeros(topsz.size - 1)
    p_value = np.zeros(topsz.size - 1)
    std_err = np.zeros(topsz.size - 1)

    for i in range(topsz.size -1 ):
        zx= z[np.where((z >=topsz[i]) & (z< topsz[i+1]))]
        vix = vi[np.where((z >=topsz[i]) & (z < topsz[i+1]))]
        slope[i], intercept[i],r_value[i],p_value[i],std_err[i] = st.linregress(zx,vix)
        print("%10.0f  %10.4f  %10.3f  %10.3f  %10.3f  %20s  %10.1f  %10.1f" %\
        (intercept[i],slope[i],r_value[i],p_value[i],std_err[i],topsname[i],min(zx),max(zx)), file=f0)
    f0.close()
    return intercept,slope,r_value,p_value,std_err

def tovel(dt, metric=False):
    """Convert dtc to vint."""
    # this is how it should be. If depth in m divide by 3.2808
    if metric:
        vi = 1000000.0 / (dt * 3.2808)
    else:
        vi = 1000000.0 / dt
    return vi

def t2d_lwz(v0,k,t1w):
    """Linear with depth velocity function."""
    depth = (v0 / k) * (np.expm1(k * t1w))
    return depth


def residuals(p,y,x):
    """."""
    err = y - peval(x,p)
    return err

def peval(t,p):
    """."""
    p[2]=0.0
    z = (t * p[1] + (t**2) * p[0])
    return z

def plot_df(zpindf,dirsplit,fname,wname,hideplot=False):
    fig,ax = plt.subplots(1,2,figsize=(8,7))
    ax[0].invert_yaxis()
    ax[0].plot(zpindf['T1W'],zpindf['Z'],c='r',label='Actual')
    ax[0].plot(zpindf['T1W'],zpindf['ZQ'],c='b',label='Quad')
    ax[0].plot(zpindf['T1W'],zpindf['ZQLSQ'],c='g',label='QuadLSQR')
    ax[0].plot(zpindf['T1W'],zpindf['ZLWZ'],c='k',label='LWZ')
    ax[0].set_xlabel('T1W in sec')
    ax[0].set_ylabel('Depth')
    ax[0].legend()
    ax[1].invert_yaxis()
    ax[1].plot(zpindf['VI'],zpindf['Z'],lw=2,c='m',label='VI')
    ax[1].plot(zpindf['VILWZ'],zpindf['Z'],lw=2,c='g',label='LWZ')
    ax[1].legend()
    ax[1].set_xlabel('VI in m/sec')
    fig.suptitle(wname)
    fig.tight_layout()

    pdfcl = os.path.join(dirsplit,fname) + ".pdf"
    if not hideplot:
        plt.show()
    fig.savefig(pdfcl)
    print(f'Successfully saved {pdfcl}')


def getcommandline():
    """."""
    parser = argparse.ArgumentParser(description='Fit Linear with Depth  and Quadratic functions to sonic or DTR. Dec 30,2018')
    parser.add_argument('datafilename',help='Depth Time 1Way file with sonic, density,shear')
    parser.add_argument('--dtcols',type=int,nargs=2,default=(0,1),help='depth  T1W columns. dfv=0 1')
    parser.add_argument('--headerlines',type=int,default=0,help='headerlines to skip. dfv=0')
    parser.add_argument('--slowness',action='store_true',default=False,
        help='Input is LAS, i.e. depth slowness. dfv= depth time pairs')
    parser.add_argument('--tometric',action='store_true',default=False,
        help='convert slowness to velocities in m/s.default= keep as ft/s')
    parser.add_argument('--nulval',type=float,default='-999.25000',help='Null value in sonic log. dfv= -999.25000')
    parser.add_argument('--dtshift',type=float, default=(0,0),nargs=2,help='depth shift time shift in ms, dfv = 0 0')
    parser.add_argument('--dtmultiplier',type=float,nargs=2,default=(1.0,1.0),help='depth and time multipliers')
    parser.add_argument('--zrange',nargs=3, type=float,default=[0,0.0,0.0],
        help='Depth range to fit functions and plot. Expected 0/1 zmin zmax values. dfv= 0 0 0')
    parser.add_argument('--qcfcorrection',nargs='+',default=[0,1,1],type=float,
        help='Input ratio to apply to quadratic coef. Used to adjust quad coef from shallow  to deep.')
    # qcfcorrection 0 0 0 means no correction to apply
    # 1 1.3 0.9 means apply correction to quadratic of least squares, only a2 and a1 is supplied, no a0
    # 2 1.3 0.9 1.5 means apply correction to quadratic with 3 coefs, a2, a1, a0
    parser.add_argument('--topsfilename',help='file name with 2 columns, depth and tops')
    parser.add_argument('--topcols',nargs=2,type=int,default=(0,1),help='columns of depth and top,dfv= 0 1')
    parser.add_argument('--topsin',choices=['time','depth'],default='depth',
        help='tops are in depth or t2w. dfv= depth')
    parser.add_argument('--topsmultshift',type=float,nargs=2,default=(1.0,0.0),
        help='Tops mutliplier and shift')
    parser.add_argument('--wellname',default='WELL-XXX',help='Well name column, dfv=WELL-XXX')
    parser.add_argument('--xycoords',nargs=2,default=(664055.00,2889245.00), type=float,
        help='X Y coordinates of well.dfv=664055.00 ,2889245.00')
    parser.add_argument('--listcoef',action='store_true',default=False,
        help='simple list of computed coefficients')
    parser.add_argument('--listdatain',action='store_true',default=False,help='List input data')
    parser.add_argument('--listdt',action='store_true',default=False,help='List d t data')
    parser.add_argument('--hideplot',action='store_true',default=False,
        help='Only save to pdf. default =show and save')

    result = parser.parse_args()
    if not result.datafilename:
        parser.print_help()
        exit()
    else:
        return result


def main():
    """Main program."""
    cmdl = getcommandline()
    if cmdl.topsfilename and cmdl.topsin == 'depth':
        ztop,tops = topsin(cmdl.topsfilename,cmdl.topcols[0],cmdl.topcols[1])
        ztop *= cmdl.topsmultshift[0]
        ztop += cmdl.topsmultshift[1]
    if cmdl.topsfilename and cmdl.topsin == 'time':
        ttop,tops = topsin(cmdl.topsfilename,cmdl.topcols[0],cmdl.topcols[1])
        ttop *= cmdl.topsmultshift[0]
        ttop += cmdl.topsmultshift[1]

    dirsplit,fextsplit = os.path.split(cmdl.datafilename)
    fname,fextn = os.path.splitext(fextsplit)
    cfname = 'allcoef.csv'

    if cmdl.slowness:
        if cmdl.dtcols[0] == 0:
            dtcolnames = ['Z','DTC']
        else:
            dtcolnames = ['DTC','Z']
        zpinsonic = pd.read_csv(cmdl.datafilename,usecols=cmdl.dtcols,header=0,names=dtcolnames,
            delim_whitespace=True,skiprows=cmdl.headerlines,na_values=float(cmdl.nulval))
        zpinsonic.dropna(inplace=True)
        zpinsonic.reset_index(drop=True, inplace=True)
        zpinsonic['DZ'] = np.ediff1d(zpinsonic['Z'],to_begin=zpinsonic['Z'][0])
        zpinsonic['DT'] = zpinsonic['DTC'] / 1000000.00
        zpinsonic['T1W'] = zpinsonic['DT'].cumsum()
        if cmdl.tometric:
            zpinsonic['VI'] = (1000000.0 / zpinsonic['DTC']) / 3.2808
        else:
            zpinsonic['VI'] = 1000000.0 / zpinsonic['DTC']

        tzcoef = np.polyfit(zpinsonic['T1W'],zpinsonic['Z'],2)

        # zi = np.polyval(tzcoef,ti)
        zpinsonic['ZQ'] = np.polyval(tzcoef,zpinsonic['T1W'])
        print('Full Depth Range: ')
        print(f'Z = {tzcoef[2]:.2f} + {tzcoef[1]:.2f} * T1W + {tzcoef[0]:.2f} * T1W^2')

        p0 = np.array(tzcoef)
        plsq = leastsq(residuals,p0,args=(zpinsonic['Z'],zpinsonic['T1W']), maxfev=2000)
        zpinsonic['ZQLSQ'] = np.polyval(plsq[0],zpinsonic['T1W'])
        print(f'Z = {plsq[0][2]:.2f} + {plsq[0][1]:.2f} * T1W + {plsq[0][0]:.2f} * T1W^2')

        zvicoef = np.polyfit(zpinsonic['Z'],zpinsonic['VI'],1)
        zpinsonic['VILWZ'] = np.polyval(zvicoef,zpinsonic['Z'])
        print(f'VINST = {zvicoef[1]:.2f} + {zvicoef[0]:.2f} * Z')
        zpinsonic['ZLWZ'] = t2d_lwz(zvicoef[1],zvicoef[0],zpinsonic['T1W'])

        if cmdl.zrange[0]:
            zpinsonicx = zpinsonic[(zpinsonic['Z'] >= cmdl.zrange[1]) & (zpinsonic['Z'] <= cmdl.zrange[2])]
            tzcoefx = np.polyfit(zpinsonicx['T1W'],zpinsonicx['Z'],2)

            # zi = np.polyval(tzcoef,ti)
            zpinsonicx['ZQ'] = np.polyval(tzcoefx,zpinsonicx['T1W'])
            print('Selected Depth Range: ')
            print(f'Z = {tzcoefx[2]:.2f} + {tzcoefx[1]:.2f} * T1W + {tzcoefx[0]:.2f} * T1W^2')

            p0x = np.array(tzcoefx)
            plsqx = leastsq(residuals,p0x,args=(zpinsonicx['Z'],zpinsonicx['T1W']), maxfev=2000)
            zpinsonicx['ZQLSQ'] = np.polyval(plsqx[0],zpinsonicx['T1W'])
            print(f'Z = {plsqx[0][2]:.2f} + {plsqx[0][1]:.2f} * T1W + {plsqx[0][0]:.2f} * T1W^2')
            tzcfr  = tzcoef / tzcoefx
            plsqr = plsq[0][:2] / plsqx[0][:2]
            print(f'\n\nQuadratic coef correction ratio: {tzcfr}')
            print(f'Quadratic coef Least Square correction ratio: {plsqr}')

            zvicoefx = np.polyfit(zpinsonicx['Z'],zpinsonicx['VI'],1)
            zpinsonicx['VILWZ'] = np.polyval(zvicoefx,zpinsonicx['Z'])
            print(f'VINST = {zvicoefx[1]:.2f} + {zvicoefx[0]:.2f} * Z')
            zpinsonic['ZLWZ'] = t2d_lwz(zvicoefx[1],zvicoefx[0],zpinsonicx['T1W'])


            pltfname = fname + f'{cmdl.zrange[1]:.0f}_{cmdl.zrange[2]:.0f}' + 'slw'
            plot_df(zpinsonicx,dirsplit,pltfname,cmdl.wellname,cmdl.hideplot)
            ztfname = fname + f'{cmdl.zrange[1]:.0f}_{cmdl.zrange[2]:.0f}_ztvix.csv'
            zpinsonicx.to_csv(ztfname,index=False)
            print(f'Successfully saved {ztfname}')
        else:
            pltfname = fname + 'td'
            plot_df(zpinsonic,dirsplit,pltfname,cmdl.wellname,cmdl.hideplot)
            ztfname = fname + 'ztvi.csv'
            zpinsonic.to_csv(ztfname,index=False)
            print(f'Successfully saved {ztfname}')

    else:
        # input data is depth one (or two) way time -> TDR
        if cmdl.dtcols[0] == 0:
            dtcolnames = ['Z','T1W']
        else:
            dtcolnames = ['T1W','Z']
        zpin = pd.read_csv(cmdl.datafilename,usecols=cmdl.dtcols,header=0,names=dtcolnames,
                delim_whitespace=True,skiprows=cmdl.headerlines,na_values=float(cmdl.nulval))
        zpin['Z'] *= cmdl.dtmultiplier[0]
        zpin['T1W'] *= cmdl.dtmultiplier[1]
        zpin['Z'] += cmdl.dtshift[0]
        zpin['T1W'] += cmdl.dtshift[1]
        zpin['DZ'] = np.ediff1d(zpin['Z'],to_begin=zpin['Z'][0])
        zpin['DT'] = np.ediff1d(zpin['T1W'],to_begin=zpin['T1W'][0])
        zpin['VI'] = zpin['DZ'] / zpin['DT']

        tzcoef = np.polyfit(zpin['T1W'],zpin['Z'],2)

        # zi = np.polyval(tzcoef,ti)
        zpin['ZQ'] = np.polyval(tzcoef,zpin['T1W'])
        print('Full Depth Range: ')
        print(f'Z = {tzcoef[2]:.2f} + {tzcoef[1]:.2f} * T1W + {tzcoef[0]:.2f} * T1W^2')

        p0 = np.array(tzcoef)
        plsq = leastsq(residuals,p0,args=(zpin['Z'],zpin['T1W']), maxfev=2000)
        zpin['ZQLSQ'] = np.polyval(plsq[0],zpin['T1W'])
        print(f'Z = {plsq[0][2]:.2f} + {plsq[0][1]:.2f} * T1W + {plsq[0][0]:.2f} * T1W^2')
        if cmdl.qcfcorrection[0] == 0:
            print('No quadratic correction to apply')
        if cmdl.qcfcorrection[0] == 1:
            plsqa = plsq[0][:2] * cmdl.qcfcorrection[1:3]
            plsq[0][:2] = plsqa
            # only a2 and a1, a0 is zero
            print('Adjusted Quadratic  Least squares Coefs:')
            print(f'Z = {plsq[0][2]:.2f} + {plsq[0][1]:.2f} * T1W + {plsq[0][0]:.2f} * T1W^2')

        if cmdl.qcfcorrection[0] == 2:
            tzcoefa = tzcoef[0][:3] * cmdl.qcfcorrection[1:4]
            tzcoef[0][:3] = tzcoefa
            # a2, a1, and a0
            print('Adjusted Quadratic Coefs:')
            print(f'Z = {tzcoef[0][2]:.2f} + {tzcoef[0][1]:.2f} * T1W + {tzcoef[0][0]:.2f} * T1W^2')

        zvicoef = np.polyfit(zpin['Z'],zpin['VI'],1)
        zpin['VILWZ'] = np.polyval(zvicoef,zpin['Z'])
        print(f'VINST = {zvicoef[1]:.2f} + {zvicoef[0]:.2f} * Z')
        # plot_df(zpin,dirsplit,fname,cmdl.hideplot)

        zpin['ZLWZ'] = t2d_lwz(zvicoef[1],zvicoef[0],zpin['T1W'])

        if cmdl.zrange[0]:
            zpinx = zpin[(zpin['Z'] >= cmdl.zrange[1]) & (zpin['Z'] <= cmdl.zrange[2])].copy()
            tzcoefx = np.polyfit(zpinx['T1W'],zpinx['Z'],2)

            # zi = np.polyval(tzcoef,ti)
            zpinx['ZQ'] = np.polyval(tzcoefx,zpinx['T1W'])
            print('Selected Depth Range: ')
            print(f'Z = {tzcoefx[2]:.2f} + {tzcoefx[1]:.2f} * T1W + {tzcoefx[0]:.2f} * T1W^2')

            p0x = np.array(tzcoefx)
            plsqx = leastsq(residuals,p0x,args=(zpinx['Z'],zpinx['T1W']), maxfev=2000)
            zpinx['ZQLSQ'] = np.polyval(plsqx[0],zpinx['T1W'])
            print(f'Z = {plsqx[0][2]:.2f} + {plsqx[0][1]:.2f} * T1W + {plsqx[0][0]:.2f} * T1W^2')
            tzcfr  = tzcoef / tzcoefx
            plsqr = plsq[0][:2] / plsqx[0][:2]
            print(f'\n\nQuadratic coef correction ratio: {tzcfr}')
            print(f'Quadratic coef Least Square correction ratio: {plsqr}')

            zvicoefx = np.polyfit(zpinx['Z'],zpinx['VI'],1)
            zpinx['VILWZ'] = np.polyval(zvicoefx,zpinx['Z'])
            print(f'VINST = {zvicoefx[1]:.2f} + {zvicoefx[0]:.2f} * Z')
            # plot_df(zpin,dirsplit,fname,cmdl.hideplot)

            zpinx['ZLWZ'] = t2d_lwz(zvicoefx[1],zvicoefx[0],zpinx['T1W'])
            pltfname = fname + f'{cmdl.zrange[1]:.0f}_{cmdl.zrange[2]:.0f}' + 'td'
            plot_df(zpinx,dirsplit,pltfname,cmdl.wellname,cmdl.hideplot)
            ztfname = fname + f'{cmdl.zrange[1]:.0f}_{cmdl.zrange[2]:.0f}_ztvix.csv'
            zpinx.to_csv(ztfname,index=False)
            print(f'Successfully saved {ztfname}')
        else:
            pltfname = fname + 'td'
            plot_df(zpin,dirsplit,pltfname,cmdl.wellname,cmdl.hideplot)
            ztfname = fname + 'ztvi.csv'
            zpin.to_csv(ztfname,index=False)
            print(f'Successfully saved {ztfname}')

    if not os.path.isfile(cfname):
        fo = open(cfname,'a+')
        headlist = ['Well','X','Y','QA0','QA1','QA2','LSQRQA0','LSQRQA1','LSQRQA2']
        print(','.join(headlist),file=fo)
        if cmdl.zrange[0]:
            cflist = [f'{cmdl.wellname}',f'{cmdl.xycoords[0]:.2f}',f'{cmdl.xycoords[1]:.2f}',f'{tzcoefx[2]:.2f}',
            f'{tzcoefx[1]:.2f}',f'{tzcoefx[0]:.2f}',f'{plsqx[0][2]:.2f}',f'{plsqx[0][1]:.2f}',f'{plsqx[0][0]:.2f}']
            print(','.join(cflist),file=fo)
        else:
            cflist = [f'{cmdl.wellname}',f'{cmdl.xycoords[0]:.2f}',f'{cmdl.xycoords[1]:.2f}',f'{tzcoef[2]:.2f}',
            f'{tzcoef[1]:.2f}',f'{tzcoef[0]:.2f}',f'{plsq[0][2]:.2f}',f'{plsq[0][1]:.2f}',f'{plsq[0][0]:.2f}']
            print(','.join(cflist),file=fo)

    else:
        fo = open(cfname,'a+')
        if cmdl.zrange[0]:
            cflist = [f'{cmdl.wellname}',f'{cmdl.xycoords[0]:.2f}',f'{cmdl.xycoords[1]:.2f}',f'{tzcoefx[2]:.2f}',
            f'{tzcoefx[1]:.2f}',f'{tzcoefx[0]:.2f}',f'{plsqx[0][2]:.2f}',f'{plsqx[0][1]:.2f}',f'{plsqx[0][0]:.2f}']
            print(','.join(cflist),file=fo)
        else:
            cflist = [f'{cmdl.wellname}',f'{cmdl.xycoords[0]:.2f}',f'{cmdl.xycoords[1]:.2f}',f'{tzcoef[2]:.2f}',
            f'{tzcoef[1]:.2f}',f'{tzcoef[0]:.2f}',f'{plsq[0][2]:.2f}',f'{plsq[0][1]:.2f}',f'{plsq[0][0]:.2f}']
            print(','.join(cflist),file=fo)

    fo.close()
    print(f'Successfully wrote to {cfname}')

if __name__ == '__main__':
    main()
