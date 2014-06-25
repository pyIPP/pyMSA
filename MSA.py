#!/usr/bin/env python

# ssh qa11
# setenv PYTHONPATH /afs/ipp/aug/ads-diags/common/python/lib
# module load python27/basic
try:
    import os, sys
    import numpy as np
    import dd_20140409 as dd, ww_20140403 as ww
    # modified mpl necessary
    sys.path.insert(1, '/afs/ipp/home/a/abock/lib/python2.7/site-packages/matplotlib-1.3.1-py2.7-linux-x86_64.egg') 
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import argparse
    from matplotlib.mlab import specgram
    from scipy.interpolate import interp1d
    import warnings
    import random
    import getpass
    from IPython import embed
    from lib.RzAmaker import makeRzAs
except Exception, e:
    print '    *** Error: %s'%e
    print '''
    This program runs best on the IPP Linux machines.
    Please execute the following commands to connect to one and
    load the necessary prerequisites:

    ssh qa01
    setenv PYTHONPATH /afs/ipp/aug/ads-diags/common/python/lib
    module load python27/basic
    ./MSA.py -h

    to try again.

    Otherwise contact abock or submit an issue at 
    https://github.com/pyIPP/pyMSA/issues'''
    sys.exit()

tmp = '''
Black
Red
Blue
Navy
Fuchsia
Purple
Maroon
Olive
Green
Lime
Aqua
Silver
Gray
Yellow
Teal
'''
mpl.rcParams['axes.color_cycle'] = tmp.split()

class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class MSAwriter(object):
    """Writes MSA shotfiles."""
        
    def __init__(self):
        super(MSAwriter, self).__init__()

    def readMSA(self, exp, shot):
        src = dd.shotfile()
        if not src.Open('MSA', shot, experiment=exp):
            return False


        # load R Z A* dR dZ dZdR from old shotfile
        R    = src.GetObject('R').data
        Z    = src.GetObject('Z').data
        A1   = src.GetObject('A1').data
        A2   = src.GetObject('A2').data
        A3   = src.GetObject('A3').data
        A4   = src.GetObject('A4').data
        A5   = src.GetObject('A5').data
        A6   = src.GetObject('A6').data
        A7   = src.GetObject('A7').data
        A8   = src.GetObject('A8').data
        A9   = src.GetObject('A9').data 
        A10  = src.GetObject('A10').data 
        dR   = src.GetObject('dR').data
        dZ   = src.GetObject('dZ').data
        dZdR = src.GetObject('dZdR').data

        # same for timebase
        TMSA = src.GetObject('T-MSA').data

        # now come g_m, g_m2, err_g_m, err_g_m2
        gm = src.GetSignalGroup('g_m')
        gm2 = src.GetSignalGroup('g_m2')
        errgm = src.GetSignalGroup('err_g_m')
        errgm2 = src.GetSignalGroup('err_g_m2')

        vbeam = src.GetParameter('misc', 'v_beam')

        pisigma = src.GetParameter('misc', 'pi/sigma')

        toReturn = (R, Z, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, dR, dZ, dZdR,
                    TMSA, gm, gm2, errgm, errgm2, vbeam, pisigma)
        # return False if something's missing
        return toReturn if None not in toReturn else False 

    def writeMSA(self, exp, shot, data):
        R, Z, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, dR, dZ, dZdR, \
        TMSA, gm, gm2, errgm, errgm2, vbeam, pisigma = data

        # todo: allow other R z A from file

        dest = ww.shotfile()
        dest.Open(experiment=exp, diagnostic='MSA', shotnumber=shot)
        #if not (dest.SetAreabase('R',  1, R) # old ww version has no feedback when writing...  
        #    and dest.SetAreabase('Z',  1, Z) # update when ww is updated
        #    ...
        #    and dest.SetParameter('misc', 'v_beam', vbeam)):
        #   dest.Close()
        #   return False
        dest.SetAreabase('R',  1, R.astype(np.float32))
        dest.SetAreabase('Z',  1, Z.astype(np.float32))
        dest.SetAreabase('A1', 1, A1.astype(np.float32))
        dest.SetAreabase('A2', 1, A2.astype(np.float32))
        dest.SetAreabase('A3', 1, A3.astype(np.float32))
        dest.SetAreabase('A4', 1, A4.astype(np.float32))
        dest.SetAreabase('A5', 1, A5.astype(np.float32))
        dest.SetAreabase('A6', 1, A6.astype(np.float32))
        dest.SetAreabase('A7', 1, A7.astype(np.float32))
        dest.SetAreabase('A8', 1, A8.astype(np.float32))
        dest.SetAreabase('A9', 1, A9.astype(np.float32))
        dest.SetAreabase('A10', 1, A10.astype(np.float32))
        dest.SetAreabase('dR', 1, dR.astype(np.float32))
        dest.SetAreabase('dZ', 1, dZ.astype(np.float32))
        dest.SetAreabase('dZdR', 1, dZdR.astype(np.float32))
        dest.SetTimebase('T-MSA', TMSA.astype(np.float32))
        dest.SetSignalGroup('g_m', gm.astype(np.float32))
        dest.SetSignalGroup('g_m2', gm2.astype(np.float32))
        dest.SetSignalGroup('err_g_m', errgm.astype(np.float32))
        dest.SetSignalGroup('err_g_m2', errgm2.astype(np.float32))
        dest.SetParameter('misc', 'v_beam', np.float32(vbeam))
        dest.SetParameter('misc', 'pi/sigma', pisigma)

        dest.Close()
        return True

    def latestMSAfile(self, shot=0, exp='AUGD'):
        return dd.PreviousShot('MSA', shot, exp)

    def latestMSCfile(self, shot=0, exp='MSED'):
        return dd.PreviousShot('MSC', shot, exp)

    def _movingAverage(self, signal, n):
        ret = np.cumsum(signal, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n

    def smooth(self, gmt, gm, smooth_window):
        dt = gmt[1] - gmt[0]
        n = max(1, int(smooth_window/dt/1e3))
        print 'smoothing over %i post-fft samples...'%n
        new_gmt = self._movingAverage(gmt, n)
        new_gm = np.array([self._movingAverage(g, n) for g in gm.T]).T
        #print 'in', gm.shape, 'out', new_gm.shape
        return new_gmt, new_gm

    pem40a, pem46a = None, None

    def _angleFromRawData(self, t, data, pem40, pem46, 
                          calib_factor=0.8, faraday_strength=0, btf=lambda x: 2.5, abs_offset=66.7,
                          nfft=2048, plot=False):
        t0, dt = t[0], t[1] - t[0]
        if self.pem40a == None or self.pem46a == None: # do only once, then cache
            self.pem40a = specgram(pem40, NFFT=nfft, Fs=1./dt, noverlap=nfft/2, mode='angle')[0]
            self.pem46a = specgram(pem46, NFFT=nfft, Fs=1./dt, noverlap=nfft/2, mode='angle')[0]

        I, f, t = specgram(data, NFFT=nfft, Fs=1./dt, noverlap=nfft/2, mode='default')
        a, f, t = specgram(data, NFFT=nfft, Fs=1./dt, noverlap=nfft/2, mode='angle')
        t += t0

        # 40-41 kHz, 46-47 kHz are the PEM harmonics we are interested in
        i40s = np.intersect1d(np.where(40e3 < f)[0], np.where(f < 41e3)[0])
        i46s = np.intersect1d(np.where(46e3 < f)[0], np.where(f < 47e3)[0])

        # indices of frequencies with maximum intensity
        i40 = i40s[np.average(I[i40s], axis=1).argmax()]
        i46 = i46s[np.average(I[i46s], axis=1).argmax()]

        # use phase difference between pem reference signal and measured signal for intensity sign information
        adiff40 = ((a[i40]-self.pem40a[i40])%(2*np.pi) > np.pi).astype(int)*2-1
        adiff46 = ((a[i46]-self.pem46a[i46])%(2*np.pi) > np.pi).astype(int)*2-1

        # use sign from previous lines to create "negative" intensities for full arctan2 usage
        I40 = np.sum(I[i40s], axis=0)**0.5 * adiff40
        I46 = np.sum(I[i46s], axis=0)**0.5 * adiff46

        if plot:
            plt.plot(t, I40)
            plt.plot(t, I46)
            plt.show()

        toReturn = -np.arctan2(calib_factor*I40, I46)*180./np.pi*0.5 # second stokes component, thus *0.5
        toReturn += faraday_strength*btf(t) + abs_offset
        toReturn[np.where(toReturn >  90)] -= 180
        toReturn[np.where(toReturn < -90)] += 180
        
        return t, toReturn

    def readMSX(self, exp, shot, nfft, use_calibration=True, upshiftPi=False):
        if use_calibration:
            # get latest calib factor (p0)
            MSC = dd.shotfile()
            if not MSC.Open('MSC', self.latestMSCfile(shot, 'MSED'), 'MSED'):
                return False
            calib_factor = -MSC.GetParameter('C_MSX', 'p0')
            # get faraday rotation degree
            faraday_degree    = MSC.GetParameter('C_Farada', 'BTF')  
            if faraday_degree == None:
                return False
            # get b1 absolute offset
            b1 = MSC.GetParameter('C_Angle', 'b1')
            MSC.Close(); del MSC
            # get BTF
            if shot < 30160:
                TOT = dd.shotfile()
                if not TOT.Open('TOT', shot, 'AUGD'):
                    return False
                BTFt = TOT.GetTimebase('BTF') 
                BTFd = -TOT.GetSignal('BTF')
                TOT.Close(); del TOT
            else:
                MBI = dd.shotfile()
                if not MBI.Open('MBI', shot, 'AUGD'):
                    return False
                BTFt = MBI.GetTimebase('BTF')
                BTFd = MBI.GetSignal('BTF')
                MBI.Close()
            BTF = interp1d(BTFt, BTFd, fill_value=np.average(BTFd), bounds_error=False)            
        else:
            calib_factor = [0.8]*10
            faraday_degree = [-0.7]*10
            BTF = interp1d([-10, 100], [-2.5, -2.5], fill_value=-2.5, bounds_error=False)
            b1 = 66.7

        src = dd.shotfile()
        if not src.Open('MSX', shot, exp):
            return False
        sg1 = src.GetSignalGroup('SG-1')
        sg2 = src.GetSignalGroup('SG-2')
        sg1t = src.GetTimebase('SG-1')/1e9 # nanoseconds...

        mse = [] # result

        # reference signals from PEMs
        pem40, pem46 = sg1[:,15], sg1[:,14]

        labels = ['pi' if l == 1 else 'sigma' if l != 0 else '' for l in src.GetParameter('CH-SETUP', 'PI/SIGMA')]
        nchan = len(labels) - labels.count('')

        # function to get data corresponding to MSA g_m entry <channel>        
        channelmap = [
            ('SG-1', 0, 'MSE  1'), # MSE box 1
            ('SG-1', 1, 'MSE  2'),
            ('SG-1', 2, 'MSE  3'),
            ('SG-1', 3, 'MSE  4'),
            ('SG-1', 4, 'MSE  5'),
            ('SG-1', 5, 'MSE  6'),
            ('SG-1', 6, 'MSE  7'),
            ('SG-1', 7, 'MSE  8'),
            ('SG-1', 8, 'MSE  9'),
            ('SG-1', 9, 'MSE 10'), # MSE box 10
            ('SG-2', 0, 'MER  1'), # MER box 1
            ('SG-2', 6, 'MER  2'), # MER box 2
        ]
        def getActualBox(channel, labels=np.array(labels)):
            return np.where(labels != '')[0][channel]
        def getData(channel, channelmap=channelmap):
            actualIndex = getActualBox(channel)
            group, index, l = channelmap[actualIndex]
            #print group, index
            return sg1[:,index] if group == 'SG-1' else sg2[:,index]

        for ch in range(nchan):
            print '%2i/%2i' %(ch+1, nchan)
            data = getData(ch)

            mset, cmse = self._angleFromRawData(sg1t, data, pem40, pem46, 
                                                calib_factor[getActualBox(ch)],
                                                faraday_degree[getActualBox(ch)],
                                                BTF,
                                                b1,
                                                nfft=2048)
            
            if upshiftPi:
                mse.append(cmse if labels[getActualBox(ch)] == 'sigma' else cmse+90)
            else:
                mse.append(cmse)

        # construct configuration object
        print 'generating geometric information from FARO (lib/mse2014.txt) and MSX/CH-SETUP...'
        rza = makeRzAs('lib/mse2014.txt')
        los = np.array([src.GetParameter('CH-SETUP', 'LOS-L%i'%i) for i in xrange(1,7)]).ravel()
        R = np.zeros(nchan); dR = np.zeros(nchan)
        Z = np.zeros(nchan); dZ = np.zeros(nchan)
        A = np.zeros((nchan,10))
        for i in xrange(nchan):
            R[i] = np.average(rza.R[np.where(los == getActualBox(i) + 1)])
            dR[i] = np.std(rza.R[np.where(los == getActualBox(i) + 1)])
            Z[i] = np.average(rza.z[np.where(los == getActualBox(i) + 1)])
            dZ[i] = np.std(rza.z[np.where(los == getActualBox(i) + 1)])
            A[i] = np.average(rza.Asigma[np.where(los == getActualBox(i) + 1)], axis=0) if labels[getActualBox(i)] == 'sigma' else \
                   np.average(rza.Api[   np.where(los == getActualBox(i) + 1)], axis=0)
        pisigma = src.GetParameter('CH-SETUP', 'PI/SIGMA')
        pisigma = pisigma[np.where(pisigma != 0)[0]]
        config = Bunch(R=R, z=Z, A=A, pisigma=pisigma, dR=dR, dZ=dZ,
            labels=['%s %s'%(channelmap[getActualBox(i)][2], 'pi' if pisigma[i] == 1 else 'sigma') for i in xrange(10)])
        print 'done'
        #embed()
        return (np.array(mset), np.array(mse).T, config)

    def getIndicesIncompatibleWithNBI(self, shot, t):
        NIS = dd.shotfile()
        if not NIS.Open('NIS', shot):
            return []
        p = NIS.GetSignalGroup('PNIQ')[:, :, 0]
        pt = NIS.GetTimebase('PNIQ')        
        # only timepoints where NBI3 > 1MW and rest at 0
        ind = np.intersect1d(np.where(p[:, 2] > 1000e3)[0], np.where(p[:, 0] == 0)[0])
        ind = np.intersect1d(ind, np.where(p[:, 1] == 0)[0])
        ind = np.intersect1d(ind, np.where(p[:, 3] == 0)[0])
        dt = pt[1]-pt[0] # use time between two NBI data points as limit
        pt = pt[ind]
        toRemove = []
        for i in range(len(t)):
            ddt = np.abs(t[i]-pt).min()
            if ddt > dt:
                toRemove.append(i)
        return toRemove

    def removeIncompatibleNBItimes(self, shot, t, gm):
        # remove items where NBI1,2,4 were on/NBI3 off
        toRemove = self.getIndicesIncompatibleWithNBI(shot, t)
        gm = np.delete(gm, toRemove, axis=0)
        t = np.delete(t, toRemove)
        return t, gm

    def write(self, args, onlyNBI3=True, nfft=2048, showPlot=False, smooth_window=None):
        res = self.readMSX(args.src_exp, args.src_num, nfft, args.use_calibration)
        if res == False:
            return False
        gmt, gm, cfg = res

        if args.channel_order != None:
            gm = gm[:, args.channel_order]

        if smooth_window != None:
            gmt, gm = self.smooth(gmt, gm, smooth_window)

        if onlyNBI3:
            gmt, gm = self.removeIncompatibleNBItimes(shot, gmt, gm)

        print 'calculated MSE data shape (time, channels):', gm.shape
        if showPlot:
            plt.plot(gmt, gm)
            plt.show()

        res = [None]*22

        # now set v_beam:
        NIS = dd.shotfile()
        if not NIS.Open('NIS', args.src_num):
            return False
        vbeam = np.sqrt(2*NIS.GetParameter('INJ1', 'UEXQ')[2]*1e3*1.602e-19 / 3.344e-27) # in m/s
        res[0] = cfg.R
        res[1] = cfg.z
        res[2] = cfg.A[:,0]
        res[3] = cfg.A[:,1]
        res[4] = cfg.A[:,2]
        res[5] = cfg.A[:,3]
        res[6] = cfg.A[:,4]
        res[7] = cfg.A[:,5]
        res[8] = cfg.A[:,6]
        res[9] = cfg.A[:,7]
        res[10] = cfg.A[:,8]
        res[11] = cfg.A[:,9]
        res[12] = cfg.dR
        res[13] = cfg.dZ
        res[14] = np.array([0.]*10)
        res[15] = gmt # TMSA
        res[16] = gm # gm
        res[17] = gm*0. # gm2
        res[18] = gm*0 + 0.2 # errgm
        res[19] = gm*0 + 180. # errgm2
        res[20] = vbeam # vbeam
        res[21] = cfg.pisigma #pisigma

        return self.writeMSA(args.dest_exp, args.dest_num if args.dest_num != None else args.src_num, res)

    def plot(self, what, args, nfft=2048, smooth_window=None):
        ax = plt.subplot(111)
        channels2use = args.only_channels
        if what=='MSA':
            res = self.readMSA(args.src_exp, args.src_num)
            piCh = np.intersect1d(np.where(res[21]==1)[0], channels2use)
            sigCh = np.intersect1d(np.where(res[21]==2)[0], channels2use)
            if res == False:
                return False
            gmt, gm = res[15:17]
            if smooth_window != None:
                gmt, gm = self.smooth(gmt, gm, smooth_window)
            lineObjects = plt.plot(gmt, gm[:,piCh], '--', dashes=(8,2)) + plt.plot(gmt, gm[:,sigCh])
            labels = ['Ch %i (pi)'%(i+1) for i in piCh] + ['Ch %i (sigma)'%(i+1) for i in sigCh]
        else:
            res = self.readMSX(args.src_exp, args.src_num, nfft, args.use_calibration, args.upshift_pi)
            if res == False:
                return False
            gmt, gm, cfg = res
            if args.channel_order != None:
                gm = gm[:, args.channel_order]
            if smooth_window != None:
                gmt, gm = self.smooth(gmt, gm, smooth_window)
            #piCh = np.intersect1d(np.where(cfg.pisigma==1)[0], channels2use)
            #sigCh = np.intersect1d(np.where(cfg.pisigma==2)[0], channels2use)
            #lineObjects = plt.plot(gmt, gm[:,piCh], '--', dashes=(8,2)) + plt.plot(gmt, gm[:,sigCh])
            #labels = ['Ch %i (pi)'%(i+1) for i in piCh] + ['Ch %i (sigma)'%(i+1) for i in sigCh]
            lineObjects = plt.plot(gmt, gm)
            labels = cfg.labels
            

        leg = plt.legend(lineObjects, [l for l in labels], fontsize=10)
        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)
        plt.title('%s:%s %i' % (args.src_exp, what, args.src_num))
        plt.xlabel('t [s]')
        plt.ylabel('g_m [deg]')
        
        from matplotlib.ticker import MultipleLocator
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.4))
        
        plt.grid(True, 'both')
        
        badind = self.getIndicesIncompatibleWithNBI(args.src_num, gmt)
        upper = np.zeros(len(gmt))
        lower = np.zeros(len(gmt))
        ax = plt.gca()
        upper[badind] = max(ax.get_ylim())
        lower[badind] = min(ax.get_ylim())
        plt.fill_between(gmt, upper, lower, color='r', alpha=0.2)

        plt.show()
        return True

    def replicate(self, args):
        print 'replicating %s:%s:%i to %s:%s:%i...' % (args.src_exp,  'MSA', args.src_num,
                                                       args.dest_exp, 'MSA', args.dest_num)
        res = self.readMSA(args.src_exp, args.src_num)
        if res == False:
            return False
        return self.writeMSA(args.dest_exp, args.dest_num, res)

def main():
    possibleActions = ['replicate','write', 'plotMSA', 'plotMSX']
    paexpl = ['re-write existing MSA file', 'write MSA file from MSX file',
        'read and plot MSA file', 'compute angles from MSX and only plot']

    parser = argparse.ArgumentParser(description='Handle MSX data.')
    parser.add_argument('-sn', '--src-num', type=int, required=True,
        help='source shotnumber, e.g. -se 29761')
    parser.add_argument('-se', '--src-exp', type=str, default='AUGD',
        help='source experiment, e.g. -se AUGD')
    parser.add_argument('-de', '--dest-exp', type=str, default=getpass.getuser().upper(),
        help='destination experiment, e.g. -de ABOCK')
    parser.add_argument('-dn', '--dest-num', type=int, default=None,
        help='destination shotnumber, e.g. -dn 123')
    parser.add_argument('-nc', '--no-calibration', dest='use_calibration', action='store_false', default=True,
        help="don't load calibration from MSC, use average hardcoded values instead")
    parser.add_argument('-i124', '--ignore-NBI124', dest='ignore_NBI124', action='store_true', default=False,
        help="don't remove data where wrong NBI configuration was present")
    parser.add_argument('-oc', '--only-channels', dest='only_channels', type=str, default=None,
        help="only plot channels 1,2,5")
    parser.add_argument('-up', '--upshift-pi', dest='upshift_pi', action='store_true', default=False,
        help="shift pi lines up by 90deg")
    #parser.add_argument('-nfft', type=int, default=2048,
    #    help='FFT window length, e.g. -nfft 2048') # can destroy coherence in phase calculation, careful when choosing other values
    nfft = 2048
    parser.add_argument('-co', '--channel-order', dest='channel_order', type=str, default='1,2,3,4,5,6,7,8,9,10',
        help='reorder channels, e.g. "-co 1,2,3,4,5,6,7,9,8,10" for shots before 30992')
    
    parser.add_argument('-s', '--smooth', dest='smooth_window', type=float, default=None,
        help='smooth result over x ms, e.g. -s 4 for 4ms moving average')
    
    parser.add_argument('action', help='one of: %s'%(
        '; '.join(['"%s" -> %s'%i for i in zip(possibleActions, paexpl)])))

    args = parser.parse_args()
    
    try:
        if args.channel_order != None:
            args.channel_order = np.array(args.channel_order.split(','), dtype=int) - 1
    except Exception, e:
        raise e

    try:
        if args.only_channels != None:
            args.only_channels = np.array(args.channel_order.split(','), dtype=int) - 1
        else:
            args.only_channels = range(10)
    except Exception, e:
        raise e


    if args.action not in possibleActions:
        print 'action must be one of these:', possibleActions
        return False

    writer = MSAwriter()

    if args.action == 'replicate':
        if not writer.replicate(args):
            print 'something went wrong :-('
            return False
    elif args.action == 'write':
        if not writer.write(args, onlyNBI3=(not args.ignore_NBI124), smooth_window=args.smooth_window):
            print 'something went wrong :-('
            return False
    elif args.action in ('plotMSA', 'plotMSX'):
        if not ((args.action == 'plotMSA' and writer.plot('MSA', args, nfft, smooth_window=args.smooth_window))
            or  (args.action == 'plotMSX' and writer.plot('MSX', args, nfft, smooth_window=args.smooth_window))):
            print 'something went wrong :-('
            return False
    return True
    

if __name__ == '__main__':
    main()


# os.system('./MSA.py -sn 29761 -se AUGD -de ABOCK -dn 234 write')
# ./MSA.py -sn 29761 -se ABOCK -de ABOCK -dn 123 replicate
# ./MSA.py -sn 29761 -se AUGD -de ABOCK -dn 234 write
