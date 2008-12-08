import healpy as H
from numpy import *
from pylab import *

nside = 32
npix = H.nside2npix(nside)

vec = array(H.pix2vec(nside,arange(npix)))
d = 3*array([0.5,0.3,0.2])
#d /= norm(d)
c = 2.

sig = dot(vec.T,d)+c
sig += randn(npix)
sig[randint(sig.size,size=10000)] = H.UNSEEN

thb,phb = H.pix2ang(nside,arange(npix))
glat=90.-thb*180/pi

w=abs(glat)<10
sig[w] += 3.*randn(w.sum())+10

H.mollview(sig)

print 'ini: mono=%s, dipole=%s'%(c,d)
mono,dipole = H.pixelfunc.fit_dipole(sig)
print 'fit: mono=%s, dipole=%s'%(mono,dipole)

resfit = dot(vec.T,dipole)+mono

diff = sig.copy()
diff[sig != H.UNSEEN] -= resfit[sig != H.UNSEEN]

H.mollview(diff)

# use remove_dipole
m2 = H.pixelfunc.remove_dipole(sig)
m3 = H.pixelfunc.remove_monopole(sig)

H.mollview(m2)
H.mollview(m3)

