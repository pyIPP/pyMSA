#!/usr/bin/env python

import os, sys
import numpy as np
import dd_local as dd, ww_local as ww
import matplotlib.pyplot as plt
import argparse
from pylab import specgram
from scipy.interpolate import interp1d

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
		Z    = src.GetObject('R').data
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
		#	dest.Close()
		#	return False
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

	def readMSX(self, exp, shot, nfft):
		# get latest calib factor (p0)
		MSC = dd.shotfile()
		if not MSC.Open('MSC', self.latestMSCfile(shot, 'MSED'), 'MSED'):
			return False
		calib_factor = -MSC.GetParameter('C_MSX', 'p0')
		# get faraday rotation degree
		faraday_degree    = MSC.GetParameter('C_Farada', 'BTF')  
		# get b1 absolute offset
		b1 = MSC.GetParameter('C_Angle', 'b1')
		MSC.Close(); del MSC
		# get BTF
		TOT = dd.shotfile()
		if not TOT.Open('TOT', shot, 'AUGD'):
			return False
		BTFt = TOT.GetTimebase('BTF') 
		BTFd = -TOT.GetSignal('BTF')
		BTF = interp1d(BTFt, BTFd, fill_value=np.average(BTFd), bounds_error=False)
		TOT.Close(); del TOT

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

	def write(self, args, onlyNBI3=True, nfft=8192, showPlot=False):

		res = self.readMSX(args.src_exp, args.src_num, nfft)
		if res == False:
			return False
		gmt, gm = res

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
		res = self.readMSA('AUGD', latestMSA)
		if res == False:
			return False
		res = list(res)

		# now set v_beam:
		t0 = gmt.min()
		t1 = gmt.max()
		nis = dd.shotfile()
		nis.Open('NIS', args.src_num)
		pniq = nis.GetSignalGroup('PNIQ')
		pniqt = nis.GetTimebase('PNIQ')
		index0 = np.abs(pniqt-t0).argmin()
		index1 = np.abs(pniqt-t1).argmin()
		vbeam = np.average(pniq[index0:index1+1, 2, 0])

		res[14] = gmt # TMSA
		res[15] = gm # gm
		res[16] = gm*0. # gm2
		res[17] = gm*0 + 0.2 # errgm
		res[18] = gm*0 + 180. # errgm2
		res[19] = vbeam # vbeam

		return self.writeMSA(args.dest_exp, args.dest_num, res, RzAfile=args.RzA_file)

	def plot(self, what, args, nfft=8192):
		if what=='MSA':
			res = self.readMSA(args.src_exp, args.src_num)
			if res == False:
				return False
			plt.plot(res[14], res[15])
			plt.show()
		else:
			res = self.readMSX(args.src_exp, args.src_num, nfft)
			if res == False:
				return False
			gmt, gm = res
			plt.plot(gmt, gm)
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
	possibleActions = ('replicate','write', 'plotMSA', 'plotMSX')

	parser = argparse.ArgumentParser(description='Write MSA shotfiles.')
	parser.add_argument('action', help='one of: %s'%(' '.join(possibleActions)))
	parser.add_argument('-se', '--src-exp', type=str, default='AUGD')
	parser.add_argument('-sn', '--src-num', type=int, required=True)
	parser.add_argument('-de', '--dest-exp', type=str, default='AUGD')
	parser.add_argument('-dn', '--dest-num', type=int)
	parser.add_argument('-rza', '--RzA-file', type=str, default=None)

	args = parser.parse_args()

	if args.action not in possibleActions:
		print 'action must be one of these:', possibleActions
		return False

	writer = MSAwriter()

	if args.action == 'replicate':
		if not writer.replicate(args):
			print 'something went wrong :-('
			return False
	elif args.action == 'write':
		if not writer.write(args):
			print 'something went wrong :-('
			return False
	elif args.action in ('plotMSA', 'plotMSX'):
		if not ((args.action == 'plotMSA' and writer.plot('MSA', args))
			or  (args.action == 'plotMSX' and writer.plot('MSX', args))):
			print 'something went wrong :-('
			return False
	print 'okey-dokey'
	return True
	

if __name__ == '__main__':
	main()


# os.system('./MSA.py -sn 29761 -se AUGD -de ABOCK -dn 234 new')
# ./MSA.py -sn 29761 -se ABOCK -de ABOCK -dn 123 replicate
# ./MSA.py -sn 29761 -se AUGD -de ABOCK -dn 234 new
