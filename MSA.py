#!/usr/bin/env python

# ssh qa11
# setenv PYTHONPATH /afs/ipp/aug/ads-diags/common/python/lib
# module load python27/basic
try:
    import os, sys
    import numpy as np
    import dd_20140409 as dd, ww_20140403 as ww
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import argparse
    from pylab import specgram
    from scipy.interpolate import interp1d
    import warnings
    import random
    import getpass
except Exception, e:
    print '(%s)'%e
    print '''This program runs best on the Linux machines.
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

# unfinished work: better plotting colors
#base = ['00', 'ff', '33']
#colors = []
#for a in base:
#	for b in base:
#		for c in base:
#			colors.append('#%s%s%s'%(a, b, c))
#random.shuffle(colors)		
#mpl.rcParams['axes.color_cycle'] = colors
#from matplotlib import cm
#sys.exit()

class MSAwriter(object):
    """Writes MSA shotfiles."""
    
    n_chan = 10
    
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

        toReturn = (R, Z, A1, A2, A3, A4, A5, A6, A7, A8, A9, dR, dZ, dZdR,
                    TMSA, gm, gm2, errgm, errgm2, vbeam)
        # return False if something's missing
        return toReturn if None not in toReturn else False 

    def writeMSA(self, exp, shot, data, RzAfile=None):
        R, Z, A1, A2, A3, A4, A5, A6, A7, A8, A9, dR, dZ, dZdR, \
        TMSA, gm, gm2, errgm, errgm2, vbeam = data

        # use new RzAs if filename given
        if RzAfile != None:
            lines = open(os.path.realpath(os.path.expanduser(RzAfile))).readlines()[1:]
            for line in lines:
                s = line.split()
                ch, nR, nZ, nAs = int(s[0])-1, float(s[1]), float(s[2]), s[3:]
                R[ch] = nR
                Z[ch] = nZ
                A1[ch], A2[ch], A3[ch], A4[ch], A5[ch], A6[ch], A7[ch], A8[ch], A9[ch] = nAs

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
        dest.SetAreabase('dR', 1, dR.astype(np.float32))
        dest.SetAreabase('dZ', 1, dZ.astype(np.float32))
        dest.SetAreabase('dZdR', 1, dZdR.astype(np.float32))
        dest.SetTimebase('T-MSA', TMSA.astype(np.float32))
        dest.SetSignalGroup('g_m', gm.astype(np.float32))
        dest.SetSignalGroup('g_m2', gm2.astype(np.float32))
        dest.SetSignalGroup('err_g_m', errgm.astype(np.float32))
        dest.SetSignalGroup('err_g_m2', errgm2.astype(np.float32))
        dest.SetParameter('misc', 'v_beam', np.float32(vbeam))
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

    def readMSX(self, exp, shot, nfft, use_calibration=True):
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
        sg1t = src.GetTimebase('SG-1')/1e9 # nanoseconds...
        mse = []
        for ch in range(self.n_chan):
            print '%2i/%2i' %(ch+1, self.n_chan)
            res = specgram(sg1[:,ch], NFFT=nfft, Fs=1./(sg1t[1]-sg1t[0]), noverlap=nfft/2)
            plt.clf() # specgram produces a plot, throw it away
            
            I = np.array(res[0])
            f = res[1]
            mset = res[2] + sg1t[0]
            
            i40 = np.intersect1d(np.where(40e3 < f)[0], np.where(f < 41e3)[0])
            i46 = np.intersect1d(np.where(46e3 < f)[0], np.where(f < 47e3)[0])

            I1 = np.sum(I[i40], axis=0)**0.5
            I2 = np.sum(I[i46], axis=0)**0.5
            #plt.plot(I1); plt.plot(I2); plt.plot(I1/I2); plt.show(); return False
            cmse = np.arctan2(calib_factor[ch]*I1, I2)*90./np.pi
            cmse += faraday_degree[ch]*BTF(mset)
            cmse += b1
            #print faraday_degree[ch], BTF(0.), b1
            mse.append(cmse-90)


        return (np.array(mset), np.array(mse).T)

    def write(self, args, onlyNBI3=True, nfft=8192, showPlot=False, smooth_window=None):

        res = self.readMSX(args.src_exp, args.src_num, nfft, args.use_calibration)
        if res == False:
            return False
        gmt, gm = res

        if args.channel_order != None:
            gm = gm[:, args.channel_order]

        if smooth_window != None:
            gmt, gm = self.smooth(gmt, gm, smooth_window)

        if onlyNBI3:
            # remove items where NBI1,2,4 were on/NBI3 off
            NIS = dd.shotfile()
            if not NIS.Open('NIS', args.src_num):
                return False
            p = NIS.GetSignalGroup('PNIQ')[:, :, 0]
            pt = NIS.GetTimebase('PNIQ')        
            # only timepoints where NBI3 > 1MW and rest at 0
            ind = np.intersect1d(np.where(p[:, 2] > 1000e3)[0], np.where(p[:, 0] == 0)[0])
            ind = np.intersect1d(ind, np.where(p[:, 1] == 0)[0])
            ind = np.intersect1d(ind, np.where(p[:, 3] == 0)[0])
            dt = pt[1]-pt[0] # use time between two NBI data points as limit
            pt = pt[ind]
            toRemove = []
            for i in range(len(gmt)):
                ddt = np.abs(gmt[i]-pt).min()
                if ddt > dt:
                    toRemove.append(i)
            gm = np.delete(gm, toRemove, axis=0)
            gmt = np.delete(gmt, toRemove)

        print 'calculated MSE data shape (time, channels):', gm.shape
        if showPlot:
            plt.plot(gmt, gm)
            plt.show()

        # get latest AUGD MSA file as a template
        latestMSA = self.latestMSAfile(args.src_num, 'AUGD')
        if latestMSA == 0:
            warnings.warn("No AUGD:MSA file found for %i, using 29000 instead."%args.src_num)
            latestMSA = 29000
        res = self.readMSA('AUGD', latestMSA)
        if res == False:
            return False
        res = list(res)

        # warning: this code was taken from mse_class, but is wrong.
        # now set v_beam:
        #t0 = gmt.min()
        #t1 = gmt.max()
        #nis = dd.shotfile()
        #nis.Open('NIS', args.src_num)
        #pniq = nis.GetSignalGroup('PNIQ')
        #pniqt = nis.GetTimebase('PNIQ')
        #index0 = np.abs(pniqt-t0).argmin()
        #index1 = np.abs(pniqt-t1).argmin()
        #vbeam = 0 #np.average(pniq[index0:index1+1, 2, 0])
        
        # self calculation instead:
        NIS = dd.shotfile()
        if not NIS.Open('NIS', args.src_num):
            return False
        vbeam = np.sqrt(2*NIS.GetParameter('INJ1', 'UEXQ')[2]*1e3*1.602e-19 / 3.344e-27) # in m/s

        res[14] = gmt # TMSA
        res[15] = gm # gm
        res[16] = gm*0. # gm2
        res[17] = gm*0 + 0.2 # errgm
        res[18] = gm*0 + 180. # errgm2
        res[19] = vbeam # vbeam

        return self.writeMSA(args.dest_exp, args.dest_num if args.dest_num != None else args.src_num, res, RzAfile=args.RzA_file)

    def plot(self, what, args, nfft=8192, smooth_window=None):
        if what=='MSA':
            res = self.readMSA(args.src_exp, args.src_num)
            if res == False:
                return False

            channels2use = args.only_channels

            gmt, gm = res[14:16]
            if smooth_window != None:
                gmt, gm = self.smooth(gmt, gm, smooth_window)
            plt.plot(gmt, gm[:, channels2use])
            plt.title('%s:%s %i' % (args.src_exp, what, args.src_num))
            plt.show()
        else:
            res = self.readMSX(args.src_exp, args.src_num, nfft, args.use_calibration)
            if res == False:
                return False
            gmt, gm = res
            if args.channel_order != None:
                gm = gm[:, args.channel_order]
            if smooth_window != None:
                gmt, gm = self.smooth(gmt, gm, smooth_window)
            
            channels2use = args.only_channels

            ax = plt.subplot(111)
            lineObjects = plt.plot(gmt, gm[:,channels2use])
            plt.legend(lineObjects, ['ch %i'%(i+1) for i in channels2use], fontsize=8)
            plt.title('%s:%s %i' % (args.src_exp, what, args.src_num))
            plt.xlabel('t [s]')
            plt.ylabel('g_m [deg]')
            
            from matplotlib.ticker import MultipleLocator
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_minor_locator(MultipleLocator(0.2))
            
            plt.grid(True, 'both')
            
            plt.show()
        return True

    def replicate(self, args):
        print 'replicating %s:%s:%i to %s:%s:%i...' % (args.src_exp,  'MSA', args.src_num,
                                                       args.dest_exp, 'MSA', args.dest_num)
        res = self.readMSA(args.src_exp, args.src_num)
        if res == False:
            return False
        return self.writeMSA(args.dest_exp, args.dest_num, res, RzAfile=args.RzA_file)

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
    parser.add_argument('-rza', '--RzA-file', type=str, default=None,
        help='optional file with new R,Z,A numbers, e.g. -rza RzAs2014.txt')
    parser.add_argument('-nc', '--no-calibration', dest='use_calibration', action='store_false', default=True,
        help="don't load calibration from MSC, use average hardcoded values instead")
    parser.add_argument('-i124', '--ignore-NBI124', dest='ignore_NBI124', action='store_true', default=False,
        help="don't remove data where wrong NBI configuration was present")
    parser.add_argument('-oc', '--only-channels', dest='only_channels', type=str, default=None,
        help="only plot channels 1,2,5")
    parser.add_argument('-nfft', type=int, default=2048,
        help='FFT window length, e.g. -nfft 2048')
    parser.add_argument('-co', '--channel-order', dest='channel_order', type=str, default='1,2,3,4,5,6,7,9,8,10',
        help='reorder channels, e.g. by default "-co 1,2,3,4,5,6,7,9,8,10"')
    
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
        if not ((args.action == 'plotMSA' and writer.plot('MSA', args, args.nfft, smooth_window=args.smooth_window))
            or  (args.action == 'plotMSX' and writer.plot('MSX', args, args.nfft, smooth_window=args.smooth_window))):
            print 'something went wrong :-('
            return False
    print 'okey-dokey'
    return True
    

if __name__ == '__main__':
    main()


# os.system('./MSA.py -sn 29761 -se AUGD -de ABOCK -dn 234 write')
# ./MSA.py -sn 29761 -se ABOCK -de ABOCK -dn 123 replicate
# ./MSA.py -sn 29761 -se AUGD -de ABOCK -dn 234 write
