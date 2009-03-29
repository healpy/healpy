"""This module provides the function query_disc and its helper functions.
Written by B. Crill.
"""

def ring_num(nside,z):
	"""Return the ring number given z of a point on the sphere.

	Input:
	 - nside: a power of 2
	 - z : a float within [-1,1]
	Return:
	 the ring number
	"""
	from numpy import sqrt
	twothird = 2.0 /3.0
	shift = 0.5
	iring = int( nside*(2.0-1.5*z) + shift)

	if (z>twothird):
		my_iring = int(nside*sqrt(3.0*(1-z)) + shift)
		iring = my_iring
		if (my_iring < 1): iring = 1
	if (z < -1*twothird):
		my_iring = int( nside* sqrt(3.0*(1.0+z)) + shift)
    		iring = my_iring
		if (my_iring < 1): iring = 1
    		iring = 4*nside - my_iring
	return iring

def ring2z(nside,iz):
	"""Return the z component from ring number.

	Input:
	 - nside: a power of 2
	 - iz: a ring number

	Return:
	 the z component of the ring
	"""
	dth1 = 1.0 / (3.0*float(nside)**2)
	dth2 = 2.0 / (3.0*float(nside))

	if (iz <= (nside-1)):
    		zring = 1.0  - float(iz)**2 * dth1
	elif (iz <=3*nside):
    		zring = float(2*nside-iz) * dth2
	else: zring = - 1.0 + float(4L*nside-iz)**2 * dth1

	return zring

def in_ring(nside,iz,phi0,dphi,nest=False):
	"""Compute the list of pixels in ring number iz in phi interval [phi0,phi0+dphi]

	Input:
	 - nside: a power of 2
	 - iz: ring number
	 - phi0: the starting longitude
	 - dphi: interval of longitude
	Keyword:
	 - nested: if True, return pixel number in nested scheme. Default: False (RING)
	Return:
	 - list of pixel numbers
	"""
	from numpy import pi,arange,concatenate,hstack,fabs,round
	from pixelfunc import nside2npix,ring2nest
	
	npix = nside2npix(nside)
	take_all = 0
	to_top = 0
	twopi = 2.0 * pi
	ncap  = 2*nside*(nside-1)
	listir = -1
	nir = 0

	phi_low = (phi0 - dphi) % twopi
	if (phi_low < 0): phi_low = phi_low + twopi
	
	phi_hi  = (phi0 + dphi) % twopi
	
	if (phi_hi < 0): phi_hi  = phi_hi  + twopi
	if (fabs(dphi-pi) < 1e-6): take_all = 1

	
	# equatorial region
	if ((iz >= nside) & (iz <= 3*nside)):
    		ir = iz - nside + 1
    		ipix1 = ncap + 4*nside*(ir-1)
    		ipix2 = ipix1 + 4*nside - 1
    		kshift = ir % 2
    		nr = nside*4
	else:
    		if (iz < nside):
        		ir = iz
        		ipix1 = 2*ir*(ir-1)
        		ipix2 = ipix1 + 4*ir - 1
		else:
        		ir = 4*nside - iz
        		ipix1 = npix - 2*ir*(ir+1)
        		ipix2 = ipix1 + 4*ir - 1
    		nr = ir*4
    		kshift = 1

	if (take_all == 1):
    		nir    = ipix2 - ipix1 + 1
		listir = arange(ipix1,ipix2+1,1)
	if (take_all == 0):
		shift = kshift * .5
    		ip_low = int(round (nr * phi_low / twopi - shift))
    		ip_hi  = int(round(nr * phi_hi  / twopi - shift))
    		ip_low = ip_low % nr     
    		ip_hi  = ip_hi  % nr      

		if (ip_low > ip_hi): to_top = 1
		
		ip_low = ip_low + ipix1
		ip_hi  = ip_hi  + ipix1
		

		if (to_top == 1):
		    	nir1 = ipix2 - ip_low + 1
    			nir2 = ip_hi - ipix1  + 1
    			nir  = nir1 + nir2
 			if ((nir1 > 0) & (nir2 > 0)):
        			#listir   = concatenate(arange(0,nir2,1)+ipix1, arange(0,nir1,1)+ip_low)
				list1 = arange(0,nir1,1)+ip_low
				list2 = arange(0,nir2,1)+ipix1
				listir = concatenate((list1,list2))
			else:
        			if (nir1 == 0) : listir = arange(0,nir2,1)+ipix1
				if (nir2 == 0) : listir = arange(0,nir1,1)+ip_low
      		else:
			nir = ip_hi - ip_low + 1
			listir = arange(0,nir,1)+ip_low  

	if (nest): listir = ring2nest(nside,listir)

	return listir

def query_disc(nside,v0,radius,nest=False,deg=True):
	"""Return the list of pixels within angle 'radius' from vector direction 'v0'

	Input:
	 - nside: a power of 2
	 - v0: the vector describing the direction of the center of the disc
	 - radius: the opening angle of the disc
	Keywords:
	 - nest: if True, pixel returned in nested scheme. Default: False (RING)
	 - deg: if False, radius angle expected in radian. Default: True (DEGREE)
	Return:
	 - list of pixel (as a numpy array)
	"""
	# assumes input in degrees
	
	from numpy import sqrt,sin,cos,pi,array,fabs,arccos,arcsin,size,empty,concatenate,arctan2,asarray
	from pixelfunc import vec2pix,nside2npix
	
	npix = nside2npix(nside)
	
	ang_conv = 1.0
	if (deg): ang_conv = pi/180.0
	
	cosang = cos(radius*ang_conv)

	# radius in radians
	radius_eff = radius * ang_conv

	v0 = asarray(v0)
	v0 /= sqrt((v0**2).sum())
	x0,y0,z0 = v0
	a = x0*x0 + y0*y0
	
	phi0 = 0.0
	if ((x0!= 0.0)|(y0!=0.0)): phi0 = arctan2(y0, x0)  

	cosphi0 = cos(phi0)

	rlat0 = arcsin(z0)
	rlat1 = rlat0 + radius_eff
	rlat2 = rlat0 - radius_eff

	if (rlat1 >= pi/2.0): zmax = 1.0 
	else                : zmax = sin(rlat1)
	
	irmin = max(ring_num(nside,zmax) - 1,1)

	if (rlat2 <= -pi/2.0): zmin = -1.0 
	else                 : zmin = sin(rlat2)
	
	irmax = min(ring_num(nside,zmin) + 1,4*nside-1)
	
	first = 1
	for iz in xrange(irmin,irmax+1):
		skip = 0
		z = ring2z(nside,iz)
		b = cosang - z*z0
		c = 1.0 - z*z
		cosdphi = b/sqrt(a*c)
		if ((x0 == 0)&(y0==0)):
			cosdphi = -1.0
			dphi = pi
		if (fabs(cosdphi) <= 1): dphi = arccos(cosdphi)
		else: 
			if (cosphi0 < cosdphi): skip = 1
			dphi = pi
		if (skip == 0):
	    		listir = in_ring(nside, iz, phi0, dphi,nest=nest)
			nir = size(listir)
			if (nir>0):
				if (first==1):
					work = listir.copy()
					first = 0
				else:
					work = concatenate((work,listir))
		

	return work
