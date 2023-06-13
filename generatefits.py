import os
import numpy as np
from astropy.io import fits

folder = "data/"




prihdr = fits.Header()
prihdr['INFO'] = "Ring-like corona intensity profile for alfa13 of Johannsen metric"
prihdr['AUTHORS'] = 'Sh. Riaz, A. Abdikamalov, D. Ayzenberg, C. Bambi'
prihdr['COMMENTS'] = "Makes a file containing tabulated intensity values."
prihdu = fits.PrimaryHDU(header=prihdr)
aux=[prihdu]


spins, dps, iscos = np.loadtxt("isco_a13_shafqat.dat", unpack=True)

fits_spin = np.unique(spins)
fits_a13 = dps.reshape(20, 30)
fits_a13.sort()

cola = fits.Column(name='a', format='1E',array=fits_spin)
colb = fits.Column(name='a13', format='30E',array=fits_a13)    
cols = fits.ColDefs([cola, colb])
tbhdu = fits.BinTableHDU.from_columns(cols)
aux.append(tbhdu)

# new_a13 = np.zeros([20, 15])

# for i in range(len(fits_spin)):
# 	indsm = np.argwhere(fits_a13[i] < 0).flatten()
# 	indsp = np.argwhere(fits_a13[i] >= 0).flatten()
# 	l1 = len(fits_a13[i][indsm])
# 	l2 = len(fits_a13[i][indsp])
# 	if l1 >= 7:
# 		new_a13[i] = np.append(fits_a13[i][indsm][l1 - 7:], fits_a13[i][indsp][:8])
# 	else:
# 		new_a13[i] = np.append(fits_a13[i][indsm], fits_a13[i][indsp][:8 + 7 - l1])
		




r_grid = np.arange(1, 25, 25./15.)

k = 0

for i in range(len(fits_spin)):
	spin = fits_spin[i]
    rh = 1.75 * (1. + np.sqrt(1.0 - spin * spin))
    hi = np.arange(34)
    h_grid = np.power(hi / 249., 2) * (500. - rh) + rh

	for dp in fits_a13[i]:
		for h in h_grid:
			Cols = []
			l = 0
			for r in r_grid:
				filename = folder + 'a%.05e.h_%.03e_r_%.03e_e_%.03e.a13_%.03e.a22_%.03e.a52_%.03e.dat' % (spin,h,r,0,dp,0,0)

				rbins, emis, intens = np.loadtxt(filename, unpack = True)
				# intens2 = intens

				# pfit1 = np.polyfit(rbins[1:5], intens[1:5], 1)
				# pfit2 = np.poly1d(pfit1)
				# # trff1.append(potrff1e(gs[i]))
				# intens2[0] = pfit2[rbins[0]]

				intens2 = np.zeros(intens.shape)

				pfit1 = np.polyfit(rbins[1:10], intens[1:10], 2)
				pfit2 = np.poly1d(pfit1)
				intens2[0] = pfit2(rbins[0])
				intens2[1:] = intens[1:]

				if not len(Cols):
					Cols.append(fits.Column(name='r_%d' % k, format='1E',array=rbins))
				Cols.append(fits.Column(name='intensity_%d_%d' % (k, l), format='1E',array=intens2))

				l += 1
				# col1 = fits.Column(name='r_%d' % k, format='1E',array=rbins)
				# col2 = fits.Column(name='intensity_%d' % k, format='1E',array=intens2)
	            # col3 = fits.Column(name='gmax', format='1E',array=fits_gmaxs)

			cols = fits.ColDefs(Cols)
			tbhdu = fits.BinTableHDU.from_columns(cols)
			aux.append(tbhdu)

			k += 1

thdulist = fits.HDUList(aux)
thdulist.writeto('Ring_corona_a13new.fits')
thdulist.close()
# thdulist = fits.HDUList(aux)





