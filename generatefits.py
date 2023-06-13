import os
import numpy as np
from astropy.io import fits

folder = "data/"




prihdr = fits.Header()
prihdr['INFO'] = "Lamppost corona intensity profile for an accretion disk of finite thickness"
prihdr['AUTHORS'] = 'A. Abdikamalov, C. Bambi'
prihdr['COMMENTS'] = "Makes a file containing tabulated intensity values."
prihdu = fits.PrimaryHDU(header=prihdr)
aux=[prihdu]


files = os.listdir("data/")
spins = [-0.998, -0.75, -0.5, -0.25, 0.0, 0.2, 0.35, 0.5, 0.6, 0.69, 0.77, 0.8373257, 0.8689509, 0.8939505, 0.91543156, 0.9346257, 0.9521743, 0.9684631, 0.9837458, 0.9982]
rmins = [8.643391, 7.9241004, 7.1784782, 6.402528, 5.6693025, 5.067473, 4.520775, 4.034323, 3.6376512, 3.2537637, 2.8874276, 2.6138475, 2.442025, 2.2888975, 2.1422594, 1.9949745, 1.8398409, 1.6649805, 1.4352659, 1.2274946]
Mdots = [0, 0.05, 0.1, 0.2, 0.3]
norms = [15279.471701493161, 15276.556011507395, 15277.704447002663, 15278.821526714228, 15279.235704605322, 15279.091634291992, 15277.42562267709, 15278.880186568247, 15278.351395483163, 15278.116645095313, 15278.361014946237, 15274.168824563556, 15277.075106446524, 15275.934389911306, 15275.612738793216, 15280.596280508325, 15281.759979442828, 15276.005117296218, 15274.77055351251, 15279.713049336997]

cola = fits.Column(name="spin", format="1E", array = spins)
cola = fits.ColDefs([cola])
tbhdu = fits.BinTableHDU.from_columns(cola)
aux.append(tbhdu)

# hdu = fits.open("rel_lp_table_v0.5b.fits")
# tb = hdu[1].data

# hdu = fits.open("/Users/askarabd/Documents/fudan/FITS/lp_thickness2.fits")

k = 0
# for spin in spins:
for ii in range(len(spins)):
	spin = spins[ii]
	norm = norms[ii]
	hmin1 = 1.7 * (1 + np.sqrt(1 - spin**2))
	hgrid1 = np.power(np.arange(40) / 39.0, 2) * (50.0 - hmin1) + hmin1
	hmin0 = 1.099 * (1 + np.sqrt(1 - spin**2))
	hgrid2 = np.power(np.arange(100) / 99.0,  2.5) * (50.0 - hmin0) + hmin0
	for mdot in Mdots:
	# for h in hgrid:
		Cols = []
		l = 0
		l2 = 0
		# for mdot in Mdots:
		for h in hgrid2:
			filename = "data/lp_%f_%f_%f.dat" % (spin, h, mdot)
			rdisk, intensity, emis_delta, inc_delta = np.loadtxt(filename, unpack=True, skiprows=1)
			if not len(Cols):
				Cols.append(fits.Column(name='r_%d' % k, format='1E',array=rdisk))
			Cols.append(fits.Column(name='intensity_%d_%d' % (k, l), format='1E',array=intensity / norm))
			Cols.append(fits.Column(name='emis_delta_%d_%d' % (k, l), format='1E',array=emis_delta))
			Cols.append(fits.Column(name='inc_delta_%d_%d' % (k, l), format='1E',array=inc_delta))

			# if filename not in files:
			# 	print(spin, h, mdot, filename, " is missing")
			# else:
			# 	try:
			# 		data = np.genfromtxt("data/"+filename, filling_values=-100, skiprows=1)
			# 		nans = np.where(data==-100)[0]
			# 		if len(nans):
			# 			print("There are nans in ", spin, h, mdot, filename)
			# 	except:
			# 		print("Problem in loading ", spin, h, mdot, filename)
			l+=1
		# tb = hdu[k+2].data
		# for h in hgrid1:
		# 	col_r = 'r_%d' % k
		# 	col_intens = 'intensity_%d_%d' % (k, l2)
		# 	col_emis = 'emis_delta_%d_%d' % (k, l2)
		# 	col_inc = 'inc_delta_%d_%d' % (k, l2)

		# 	if np.max(np.abs(tb[col_r] - rdisk)) > 10e-3:
		# 		print("Mismatch in ", spin, h, mdot)
		# 	intensity = tb[col_intens]
		# 	emis_delta = tb[col_emis]
		# 	inc_delta = tb[col_inc]
		# 	Cols.append(fits.Column(name='intensity_%d_%d' % (k, l), format='1E',array=intensity / norm))
		# 	Cols.append(fits.Column(name='emis_delta_%d_%d' % (k, l), format='1E',array=emis_delta))
		# 	Cols.append(fits.Column(name='inc_delta_%d_%d' % (k, l), format='1E',array=inc_delta))

		# 	l+=1
		# 	l2+=1



		cols = fits.ColDefs(Cols)
		tbhdu = fits.BinTableHDU.from_columns(cols)
		aux.append(tbhdu)
		k+=1


thdulist = fits.HDUList(aux)
thdulist.writeto('lp_thickness_v1.1b.fits')
thdulist.close()



# r_grid = np.arange(1, 25, 25./15.)

# k = 0

# for i in range(len(fits_spin)):
# 	spin = fits_spin[i]
#     rh = 1.75 * (1. + np.sqrt(1.0 - spin * spin))
#     hi = np.arange(34)
#     h_grid = np.power(hi / 249., 2) * (500. - rh) + rh

# 	for dp in fits_a13[i]:
# 		for h in h_grid:
# 			Cols = []
# 			l = 0
# 			for r in r_grid:
# 				filename = folder + 'a%.05e.h_%.03e_r_%.03e_e_%.03e.a13_%.03e.a22_%.03e.a52_%.03e.dat' % (spin,h,r,0,dp,0,0)

# 				rbins, emis, intens = np.loadtxt(filename, unpack = True)
# 				# intens2 = intens

# 				# pfit1 = np.polyfit(rbins[1:5], intens[1:5], 1)
# 				# pfit2 = np.poly1d(pfit1)
# 				# # trff1.append(potrff1e(gs[i]))
# 				# intens2[0] = pfit2[rbins[0]]

# 				intens2 = np.zeros(intens.shape)

# 				pfit1 = np.polyfit(rbins[1:10], intens[1:10], 2)
# 				pfit2 = np.poly1d(pfit1)
# 				intens2[0] = pfit2(rbins[0])
# 				intens2[1:] = intens[1:]

# 				if not len(Cols):
# 					Cols.append(fits.Column(name='r_%d' % k, format='1E',array=rbins))
# 				Cols.append(fits.Column(name='intensity_%d_%d' % (k, l), format='1E',array=intens2))

# 				l += 1
# 				# col1 = fits.Column(name='r_%d' % k, format='1E',array=rbins)
# 				# col2 = fits.Column(name='intensity_%d' % k, format='1E',array=intens2)
# 	            # col3 = fits.Column(name='gmax', format='1E',array=fits_gmaxs)

# 			cols = fits.ColDefs(Cols)
# 			tbhdu = fits.BinTableHDU.from_columns(cols)
# 			aux.append(tbhdu)

# 			k += 1

# thdulist = fits.HDUList(aux)
# thdulist.writeto('Ring_corona_a13new.fits')
# thdulist.close()
# # thdulist = fits.HDUList(aux)





