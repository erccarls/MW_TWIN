'''
Created on Apr 22, 2012

@author: carlson
'''
import os, sys
sys.path.append('./pyfits/lib')

import pyfits, numpy, math #@UnresolvedImport

from sys import stdout

file1   = sys.argv[1]
file2   = sys.argv[2]
fileout = sys.argv[3]

hdulist1 = pyfits.open(file1, mode='update')
scidata1 = hdulist1[0].data
hdulist1.info()

hdulist2 = pyfits.open(file2, mode='update')
scidata2 = hdulist2[0].data
hdulist2.info()
    
# Read FITS header to get dimensions of emissivity map (pixels)
dimen = scidata1.shape
iRange = int(dimen[0])
jRange = int(dimen[1])
kRange = int(dimen[2])
lRange = int(dimen[3])
print "Ranges are", iRange, jRange, kRange, lRange

image = numpy.zeros((iRange,jRange, kRange, lRange))
print 'Subtracting'

for j in range(0,jRange):                 # x loop
    for k in range(0,kRange):             # y loop
        for l in range(0,lRange):        # E-loop
            diff = scidata1[0,j,k, l]-scidata2[0,j,k, l]
            image[0,j,k,l] = diff
            
print 'Complete.  FITS image output to ' + fileout

#===========================================================================
# Write array to FITS image
#===========================================================================    
hdu     = pyfits.PrimaryHDU(image)
hdulist = pyfits.HDUList([hdu])
if os.access(fileout, os.F_OK ):  os.remove(fileout)
hdulist.writeto(fileout)    

