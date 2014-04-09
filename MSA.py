#!/usr/bin/env python

import os, sys
import numpy as np
import dd_local as dd, ww_local as ww
import matplotlib.pyplot as plt
import argparse
from pylab import specgram

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
		t0 = TMSA.min()
		t1 = TMSA.max()

		# now come g_m, g_m2, err_g_m, err_g_m2
		gm = src.GetSignalGroup('g_m')
		gm2 = src.GetSignalGroup('g_m2')
		errgm = src.GetSignalGroup('err_g_m')
		errgm2 = src.GetSignalGroup('err_g_m2')

		# now set v_beam:
		nis = dd.shotfile()
		nis.Open('NIS', shot)
		pniq = nis.GetSignalGroup('PNIQ')
		pniqt = nis.GetTimebase('PNIQ')
		index0 = np.abs(pniqt-t0).argmin()
		index1 = np.abs(pniqt-t1).argmin()
		vbeam = np.average(pniq[index0:index1+1, 2, 0])

		toReturn = (R, Z, A1, A2, A3, A4, A5, A6, A7, A8, A9, dR, dZ, dZdR,
			        TMSA, gm, gm2, errgm, errgm2, vbeam)
		# return False if something's missing
		return toReturn if None not in toReturn else False 

	def writeMSA(self, exp, shot, data):
		R, Z, A1, A2, A3, A4, A5, A6, A7, A8, A9, dR, dZ, dZdR, \
		TMSA, gm, gm2, errgm, errgm2, vbeam = data
		dest = ww.shotfile()
		dest.Open(experiment=exp, diagnostic='MSA', shotnumber=shot)
		#if not (dest.SetAreabase('R',  1, R) # old ww version has no feedback when writing...  
		#    and dest.SetAreabase('Z',  1, Z) # update when ww is updated
		#    ...
		#    and dest.SetParameter('misc', 'v_beam', vbeam)):
		#	dest.Close()
		#	return False
		dest.SetAreabase('R',  1, R)
		dest.SetAreabase('Z',  1, Z)
		dest.SetAreabase('A1', 1, A1)
		dest.SetAreabase('A2', 1, A2)
		dest.SetAreabase('A3', 1, A3)
		dest.SetAreabase('A4', 1, A4)
		dest.SetAreabase('A5', 1, A5)
		dest.SetAreabase('A6', 1, A6)
		dest.SetAreabase('A7', 1, A7)
		dest.SetAreabase('A8', 1, A8)
		dest.SetAreabase('A9', 1, A9)
		dest.SetAreabase('dR', 1, dR)
		dest.SetAreabase('dZ', 1, dZ)
		dest.SetAreabase('dZdR', 1, dZdR)
		dest.SetTimebase('T-MSA', TMSA)
		dest.SetSignalGroup('g_m', gm)
		dest.SetSignalGroup('g_m2', gm2)
		dest.SetSignalGroup('err_g_m', errgm)
		dest.SetSignalGroup('err_g_m2', errgm2)
		dest.SetParameter('misc', 'v_beam', np.float32(vbeam))
		dest.Close()
		return True

	def readMSX(self, exp, shot, nfft=2048):
		calib_factor = 0.8 # is that really 0.8 for all channels?
		faraday_degree    = -0.7  # where is that from?
		BTF               =  2.5  # need to get from shotfile system

		src = dd.shotfile()
		if not src.Open('MSX', shot, exp):
			return False
		sg1 = src.GetSignalGroup('SG-1')
		sg1t = src.GetTimebase('SG-1')/1e9
		mse = []
		for ch in range(self.n_chan):
			print '%2i/%2i' %(ch+1, self.n_chan)
			res = specgram(sg1[:,ch], NFFT=nfft, Fs=1./(sg1t[1]-sg1t[0]), noverlap=nfft/4)
			plt.clf() # specgram produces a plot, throw it away
			
			I = np.array(res[0])
			f = res[1]
			mset = res[2] + sg1t[0]
			
			i40 = np.intersect1d(np.where(40e3 < f)[0], np.where(f < 41e3)[0])
			i46 = np.intersect1d(np.where(46e3 < f)[0], np.where(f < 47e3)[0])

			I1 = np.sqrt(np.sum(I[i40], axis=0))
			I2 = np.sqrt(np.sum(I[i46], axis=0))

			mse.append(np.arctan2(calib_factor*I1, I2)*90./np.pi+faraday_degree*BTF)

		return (mset, np.array(mse).T)

	def new(self, args):
		res = self.readMSX(args.src_exp, args.src_num)
		if res == False:
			return False
		gmt, gm = res
		plt.plot(gmt, gm)
		plt.show()
		return True

	def replicate(self, args):
		print 'replicating %s:%s:%i to %s:%s:%i...' % (args.src_exp, args.src_diag, args.src_num,
			                                           args.dest_exp,    'MSA',     args.dest_num)
		res = self.readMSA(args.src_exp, args.src_num)
		if res == False:
			return False
		return self.writeMSA(args.dest_exp, args.dest_num, res)


def main():
	parser = argparse.ArgumentParser(description='Write MSA shotfiles.')
	parser.add_argument('action')
	parser.add_argument('-se', '--src-exp', type=str, default='AUGD')
	parser.add_argument('-sn', '--src-num', type=int, required=True)
	parser.add_argument('-de', '--dest-exp', type=str, required=True)
	parser.add_argument('-dn', '--dest-num', type=int, required=True)

	args = parser.parse_args()
	
	possibleActions = ('replicate','new')
	if args.action not in possibleActions:
		print 'action must be one of these:', possibleActions
		return False

	writer = MSAwriter()

	if args.action == 'replicate':
		if not writer.replicate(args):
			print 'something went wrong :-('
			return False
	elif args.action == 'new':
		if not writer.new(args):
			print 'something went wrong :-('
			return False
	print 'okey-dokey'
	return True
	

if __name__ == '__main__':
	main()


# os.system('./MSA.py -sn 29761 -se AUGD -de ABOCK -dn 234 new')
# ./MSA.py -sn 29761 -se ABOCK -de ABOCK -dn 123 replicate
# ./MSA.py -sn 29761 -se AUGD -de ABOCK -dn 234 new