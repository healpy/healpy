import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

SIZE = 400
DPI = 60

m = np.arange(hp.nside2npix(32))
hp.mollview(m, nest = True, xsize=SIZE)
plt.savefig('static/moll_nside32_nest.png', title='Mollview image NESTED', dpi=DPI)

hp.mollview(m, nest = False, xsize=SIZE)
plt.savefig('static/moll_nside32_ring.png', title='Mollview image RING', dpi=DPI)
