#===============================================================================
# Author: Eric Carlson 
# E-mail: erccarls@ucsc.edu
# Description: Integrates galprop emissivity maps for arbitrary inclination through LOS.
# Output:Outputs FITS with units mJy/beam for a specified frequency.S
# Last Modified: May 23 2012
#===============================================================================
from scipy import *
import sys
sys.path.append('./pyfits/lib')


#===============================================================================
# Parameters: Emissivity filename, FITS ouput filename, physical radius of
# emissivity file, and total physical height (z) of emissivity file (i.e. +-5kpc
# would be take arg 10),inclination in degrees,  physical distance to the object
# in kpc, frequency, and the instrument half power beam width (hpbw) in arcsec
# and the faintest contour n, and the width in arc-minutes of the output plot.
#===============================================================================
def calc_emiss(emissFilename, fileout, r_kpc, z_kpc, inclination, objectDist, frequency, hpbw,contour_n,width):
    import pyfits, numpy, math #@UnresolvedImport
    import os, sys
    from sys import stdout
    
    pc2cm = 3.08568025e18
    kpc2cm = 3.08568025e21
    sqArcSecPerSr  = 42545949625.0  # Square arcseconds per sr
    
    # Determine which energy bin to use based on the frequency passed in.
    vSyncLow     = 1.02026e8        #These are found in the GALDEF File
    vSyncFactor  = 1.125            #These are found in the GALDEF File
    
    # Determine correct energy bin
    eBin = int(round(math.log(frequency/vSyncLow)/math.log(vSyncFactor))) 
    print 'Selected Frequency:' + str(frequency/1e9) + 'GHz \nEnergy bin: ' + str(eBin)

    
    # generate a square image.
    outRes = 300         # Number of pixels for square output image was 300 for paper
    num_zSteps = 100     # Number of z-steps along the integration line

    # load emissivity data
    hdulist = pyfits.open(emissFilename, mode='update')
    scidata = hdulist[0].data
    hdulist.info()
    
    # Read FITS header to get dimensions of emissivity map (pixels)
    dimen = scidata.shape
    iRange = int(dimen[0])
    jRange = int(dimen[1])
    kRange = int(dimen[2])
    lRange = int(dimen[3])
    print "Ranges are", iRange, jRange, kRange, lRange

    if (eBin>=jRange): 
        print "Energy out of range.  Exiting..."
        return
    
    #===========================================================================
    # We now will choose a coordinate system of (x',y',z') with the +z' face of a
    # cube parallel to the observer (output image) so that we integrate along z'
    # axis. Thus we have a cylindrical input map rotated in a cube. This cube is
    # chosen to have side length of 1.5*the diameter of the cylinder.  Points in
    # the primed system are mapped to points in the unprimed system via an inverse
    # euler rotation.  Rotations are unitary so the transpose of the standard
    # rotation matrix yields the inverse.
    #===========================================================================
    
    # Calculate the physical side length of the output map based on object distance and angular size of output. 
    sideLength = 40         # in kpc
    print 'Output sidelength (kpc): ' , sideLength
    
    kpcPerDeltaZ   = float(sideLength)/num_zSteps            # Physical distance for z integration interval
    kpcPerPixelInZ = float(z_kpc)/kRange                     # Input Z has roughly 10:1 aspect ratio z:r
    kpcPerPixelInR = float(r_kpc)/lRange                     # 
    
    kpcPerPixelOut = float(sideLength)/float(outRes)         # Output has square aspect ratio
    cmPerPixelOut  = kpcPerPixelOut*kpc2cm                   # Output cm/pixel
    solidAngle     = 1/((objectDist*kpc2cm)**2.0)              # 4 pi * ratio of 1 cm^2 to the entire sphere given distance
      
    volPerPixelOut = (cmPerPixelOut**2.0)*kpcPerDeltaZ*kpc2cm  # one integration volume in cm^3
    
     
    srPerPixel = (kpcPerPixelOut/(objectDist))**2   # 4 pi * ratio of areas
    SqArcSecPerPixel = sqArcSecPerSr*srPerPixel     # ArcSec^2/Pixel
    pixelPerArcSecSq  = 1/(SqArcSecPerPixel)   # pixel/ArcSec^2  (Ratio of areas times 4 pi * conversion)
    beamConversion    = pixelPerArcSecSq * (math.pi * (hpbw/2)**2) / math.log(2)         # pixel/ArcSec^2  * arcSec^2/beam        
    
    #pixPerBeam =  (math.pi * (hpbw/2)**2)/SqArcSecPerPixel / math.log(2)
    #print 'ppb', pixPerBeam
    
    
    print 'ArcSec^2/pixel: ' + str(1/pixelPerArcSecSq)
    print 'Beam Conversion: ' + str(beamConversion)
    
    image = numpy.zeros((outRes, outRes))
    print 'Beginning Integration'
    
    for i in range(0,outRes):                  # x loop
        for j in range(0,outRes):              # y loop
            flux = 0.0                         # total flux sum
            for k in range (0,num_zSteps):     # z-integration
                # Shift coordinates so they are centered on galaxy and convert to kpc before transforming to emiss. coordinates.
                outPosition = ((i-outRes/2.0)*kpcPerPixelOut,(j-outRes/2.0)*kpcPerPixelOut,(k-num_zSteps/2.0)*kpcPerDeltaZ)
                inPosition  = rotate_coordinates(outPosition,inclination[0]/180.0*math.pi,0,inclination[1]/180.0*math.pi)       # Rotate to input emissivity file coords.
                r = math.sqrt(inPosition[0]**2.0+inPosition[1]**2.0)                          # Calc r in kpc
                z = inPosition[2]                                                             # Calc z in kpc
                if (abs(z) <= z_kpc/2.0 and r <= float(r_kpc)):                                 # Check that we are within the bounds of the input image.
                    xPixel = int(r/kpcPerPixelInR)                                          # Calc the actual pixel coordinates.
                    yPixel = int((z+z_kpc/2.0)/kpcPerPixelInZ)                              # Calc the actual pixel coordinates.
                    
                    # Ensure we are have a valid coordinate.  The edge flux may be underestimated since we are truncating, but assuming the emiss file is very small at boundaries this won't be an issue.
                    if (xPixel < lRange and yPixel < kRange):    
                        density = scidata[0,eBin,yPixel, xPixel]                              # Read density from file                        
                        flux += volPerPixelOut*density                                        # Volume*energy density
            if(flux!=0):            
                image[j,i] = flux*frequency*solidAngle # Image in erg/s/cm^2
        
        # Monitor Progress in terminal.
        sys.stdout.write('\r' + str(int(float(i)/outRes*100+1.0)))
        sys.stdout.flush()
    
    print 'Integration Complete...\n'
    print 'Total Luminosity (erg/s):' + str(numpy.sum(image)/solidAngle*4*math.pi)
    print 'Total Flux (Jy):' + str(numpy.sum(image)*1.0e23/frequency) # in Jy
    print 'Complete.  FITS image output to ' + fileout

    #===========================================================================
    # Write array to FITS image
    #===========================================================================    
    hdu     = pyfits.PrimaryHDU(beamConversion*image*1.0e23*1000/frequency) # in mJy/beam
    hdulist = pyfits.HDUList([hdu])
    if os.access(fileout, os.F_OK ):  os.remove(fileout)
    hdulist.writeto(fileout)    
    

    #===========================================================================
    # Contour plot with scipy
    #===========================================================================
    import matplotlib.pyplot as plt  #@UnresolvedImport
    import matplotlib.image as mpimg #@UnresolvedImport
    import matplotlib.cm as cm
    from matplotlib.patches import Ellipse
    import matplotlib.patches as mpatches
    import matplotlib
    from matplotlib.collections import PatchCollection
    
    SideLengthArcMin = math.atan(sideLength/float(objectDist))*180/math.pi*60 #side length in arc minutes

    img = plt.imshow(beamConversion*image*1.0e23/frequency*1000,origin='lower', extent=[-SideLengthArcMin/2,SideLengthArcMin/2,-SideLengthArcMin/2,SideLengthArcMin/2],vmin=0, vmax=15.0)
    #img.set_clim=(0.0,15.0)
    
    tickList = [2**-1.5,2**0,2**.5,2**1.0,2**1.5,2**2,2**2.5,2**3,2**3.5,15.0]
    for i in range(len(tickList)):
        tickList[i] = round(tickList[i],2)
    plt.colorbar(ticks=tickList, orientation='vertical')
    
    
    pltSideLengthKpc = 20.0 # kpc
    
    # Radio contour is so small we prefer a different sidelength.
    if  (fileout == 'm31_out_DM_ONLY.fits'):
        print 'M31, Adjusting sidelength'
        pltSideLengthKpc = 10.0 # kpc
    
    pltSideLength = math.atan(pltSideLengthKpc/float(objectDist))*180/math.pi*60
    plt.xlim((-pltSideLength/2,pltSideLength/2))
    plt.ylim((-pltSideLength/2,pltSideLength/2))
    
    #plt.colorbar( img, shrink=.8 ,extend='both')
    
    plt.xlabel('x (arcmin)')
    plt.ylabel('y (arcmin)')
    plt.title('Flux ($mJy Beam^{-1}$ at ' + str(frequency/1.0e9) + 'GHz)\nTotal Flux:' + str(round(numpy.sum(image)*1.0e23/frequency,6)) + ' Jy')
    levelList =  []
    for i in range(0,20):
        levelList.append(2.0**((contour_n+i)/2))
    
    CS = plt.contour(beamConversion*image*1.0e23/frequency*1000,origin='lower', extent=[-SideLengthArcMin/2,SideLengthArcMin/2,-SideLengthArcMin/2,SideLengthArcMin/2],levels=levelList, colors = 'w')
    
    plt.clabel(CS, fontsize=8, inline=1)
    #plt.colorbar( shrink=.8 ,extend='both')
    #plt.flag()

    plt.savefig(str(fileout)+ '.png')
    #plt.show()
    plt.clf()



# Rotates coordinates given full set of euler angles.  1 Redundant angle for axially symmetric galaxy
#===============================================================================
# Maps coordinates in the output image to those of the input map via inverse
# euler rotation matrix.  See wikipedia article on "Rotation Matrix".  Input is
# the rotated coordinate vector (x',y',z') followed by the euler rotation angles
# for the galaxy. Returns (x,y,z) for the emissivity map (x^2+y^2=r).
#===============================================================================
def rotate_coordinates(vector,phi,theta,psi):
    import math
    x = vector[0]*math.cos(theta)*math.cos(psi) + vector[1]*math.cos(theta)*math.sin(psi) - vector[2]*math.sin(theta)
    y = vector[0]*(-math.cos(phi)*math.sin(psi)+math.sin(phi)*math.sin(theta)*math.cos(psi)) + vector[1]*(math.cos(phi)*math.cos(psi)+math.sin(phi)*math.sin(theta)*math.sin(psi)) + vector[2]*math.sin(phi)*math.cos(theta)
    z = vector[0]*(math.sin(phi)*math.sin(psi)+math.cos(phi)*math.sin(theta)*math.cos(psi))+vector[1]*(-math.sin(phi)*math.cos(psi)+math.cos(phi)*math.sin(theta)*math.sin(psi))+vector[2]*math.cos(phi)*math.cos(theta)
    return (x,y,z)
    
    
from scipy import *
import sys
sys.path.append('./pyfits/lib')





def calcz_flux(emissFilename):
    import pyfits, numpy, math #@UnresolvedImport
    import os, sys
    
    hdulist = pyfits.open(emissFilename, mode='update')
    scidata = hdulist[0].data
    
    eBin = 23 # 1.49 ghz bin
    
    dimen = scidata.shape
    iRange = int(dimen[0])
    jRange = int(dimen[1])
    kRange = int(dimen[2])
    lRange = int(dimen[3])



#===============================================================================
# Used for contour flux comparison
#===============================================================================
def calc_emiss_limited(emissFilename, fileout, r_kpc, z_kpc, inclination, objectDist, frequency, hpbw,contour_n,width):
    import pyfits, numpy, math #@UnresolvedImport
    import os, sys
    from sys import stdout
    
    pc2cm = 3.08568025e18
    kpc2cm = 3.08568025e21
    sqArcSecPerSr  = 42545949625.0  # Square arcseconds per sr
    
    # Determine which energy bin to use based on the frequency passed in.
    vSyncLow     = 1.02026e8        #These are found in the GALDEF File
    vSyncFactor  = 1.125            #These are found in the GALDEF File
    
    # Determine correct energy bin
    eBin = int(round(math.log(frequency/vSyncLow)/math.log(vSyncFactor))) 
    print 'Selected Frequency:' + str(frequency/1e9) + 'GHz \nEnergy bin: ' + str(eBin)

    
    # generate a square image.
    outRes = 300         # Number of pixels for square output image was 300 for paper
    num_zSteps = 100     # Number of z-steps along the integration line

    # load emissivity data
    hdulist = pyfits.open(emissFilename, mode='update')
    scidata = hdulist[0].data
    hdulist.info()
    
    # Read FITS header to get dimensions of emissivity map (pixels)
    dimen = scidata.shape
    iRange = int(dimen[0])
    jRange = int(dimen[1])
    kRange = int(dimen[2])
    lRange = int(dimen[3])
    print "Ranges are", iRange, jRange, kRange, lRange

    if (eBin>=jRange): 
        print "Energy out of range.  Exiting..."
        return
    
    #===========================================================================
    # We now will choose a coordinate system of (x',y',z') with the +z' face of a
    # cube parallel to the observer (output image) so that we integrate along z'
    # axis. Thus we have a cylindrical input map rotated in a cube. This cube is
    # chosen to have side length of 1.5*the diameter of the cylinder.  Points in
    # the primed system are mapped to points in the unprimed system via an inverse
    # euler rotation.  Rotations are unitary so the transpose of the standard
    # rotation matrix yields the inverse.
    #===========================================================================
    
    # Calculate the physical side length of the output map based on object distance and angular size of output. 
    sideLength = 40         # in kpc
    print 'Output sidelength (kpc): ' , sideLength
    
    kpcPerDeltaZ   = float(sideLength)/num_zSteps            # Physical distance for z integration interval
    kpcPerPixelInZ = float(z_kpc)/kRange                     # Input Z has roughly 10:1 aspect ratio z:r
    kpcPerPixelInR = float(r_kpc)/lRange                     # 
    
    kpcPerPixelOut = float(sideLength)/float(outRes)         # Output has square aspect ratio
    cmPerPixelOut  = kpcPerPixelOut*kpc2cm                   # Output cm/pixel
    solidAngle     = 1/((objectDist*kpc2cm)**2.0)              # 4 pi * ratio of 1 cm^2 to the entire sphere given distance
      
    volPerPixelOut = (cmPerPixelOut**2.0)*kpcPerDeltaZ*kpc2cm  # one integration volume in cm^3
    
     
    srPerPixel = (kpcPerPixelOut/(objectDist))**2   # 4 pi * ratio of areas
    SqArcSecPerPixel = sqArcSecPerSr*srPerPixel     # ArcSec^2/Pixel
    pixelPerArcSecSq  = 1/(SqArcSecPerPixel)   # pixel/ArcSec^2  (Ratio of areas times 4 pi * conversion)
    beamConversion    = pixelPerArcSecSq * (math.pi * (hpbw/2)**2) / math.log(2)         # pixel/ArcSec^2  * arcSec^2/beam        
    
    #pixPerBeam =  (math.pi * (hpbw/2)**2)/SqArcSecPerPixel / math.log(2)
    #print 'ppb', pixPerBeam
    
    
    print 'ArcSec^2/pixel: ' + str(1/pixelPerArcSecSq)
    print 'Beam Conversion: ' + str(beamConversion)
    
    image = numpy.zeros((outRes, outRes))
    print 'Beginning Integration'
    
    for i in range(0,outRes):                  # x loop
        for j in range(0,outRes):              # y loop
        flux = 0.0                         # total flux sum
        for k in range (0,num_zSteps):     # z-integration
            # Shift coordinates so they are centered on galaxy and convert to kpc before transforming to emiss. coordinates.
            outPosition = ((i-outRes/2.0)*kpcPerPixelOut,(j-outRes/2.0)*kpcPerPixelOut,(k-num_zSteps/2.0)*kpcPerDeltaZ)
            inPosition  = rotate_coordinates(outPosition,inclination[0]/180.0*math.pi,0,inclination[1]/180.0*math.pi)       # Rotate to input emissivity file coords.
            r = math.sqrt(inPosition[0]**2.0+inPosition[1]**2.0)                          # Calc r in kpc
            z = inPosition[2]                                                             # Calc z in kpc
            if (abs(z) <= z_kpc/2.0 and r <= float(r_kpc)):                                 # Check that we are within the bounds of the input image.
                xPixel = int(r/kpcPerPixelInR)                                          # Calc the actual pixel coordinates.
                yPixel = int((z+z_kpc/2.0)/kpcPerPixelInZ)                              # Calc the actual pixel coordinates.
                
                # Ensure we are have a valid coordinate.  The edge flux may be underestimated since we are truncating, but assuming the emiss file is very small at boundaries this won't be an issue.
                if (xPixel < lRange and yPixel < kRange):    
                    density = scidata[0,eBin,yPixel, xPixel]                              # Read density from file                        
                    flux += volPerPixelOut*density                                        # Volume*energy density
        if(flux!=0):            
            image[j,i] = flux*frequency*solidAngle # Image in erg/s/cm^2
        
        # Monitor Progress in terminal.
        sys.stdout.write('\r' + str(int(float(i)/outRes*100+1.0)))
        sys.stdout.flush()
    
    print 'Integration Complete...\n'
    print 'Total Luminosity (erg/s):' + str(numpy.sum(image)/solidAngle*4*math.pi)
    print 'Total Flux (Jy):' + str(numpy.sum(image)*1.0e23/frequency) # in Jy
    print 'Complete.  FITS image output to ' + fileout

    #===========================================================================
    # Write array to FITS image
    #===========================================================================    
    hdu     = pyfits.PrimaryHDU(beamConversion*image*1.0e23*1000/frequency) # in mJy/beam
    hdulist = pyfits.HDUList([hdu])
    if os.access(fileout, os.F_OK ):  os.remove(fileout)
    hdulist.writeto(fileout)    






#===============================================================================
# Parameters: Emissivity filename, FITS ouput filename, physical radius of
# emissivity file, and total physical height (z) of emissivity file (i.e. +-5kpc
# would be take arg 10),inclination in degrees,  physical distance to the object
# in kpc, frequency, and the instrument half power beam width (hpbw) in arcsec
# and the faintest contour n, and the width in arc-minutes of the output plot.
#===============================================================================
def intEllipse(FITSFile, semiMaj,semiMin,dist,posAng,hpbw):
    import pyfits, numpy, math #@UnresolvedImport
    import os, sys
    from sys import stdout
    
    pc2cm = 3.08568025e18
    kpc2cm = 3.08568025e21
    sqArcSecPerSr  = 42545949625.0  # Square arcseconds per sr
    
    # Determine correct energy bin

    # load emissivity data
    hdulist = pyfits.open(FITSFile, mode='update')
    scidata = hdulist[0].data
    #hdulist.info()
    
    # Read FITS header to get dimensions of emissivity map (pixels)
    dimen = scidata.shape
    
    jRange = int(dimen[0])
    kRange = int(dimen[1])
    #print "Ranges are", jRange, kRange
    
    sideLength = 40.0         # in kpc
    SideLengthArcMin = math.atan(sideLength/float(dist))*180./math.pi*60. #side length in arc minutes

    # Calculate the physical side length of the output map based on object distance and angular size of output. 
    
    #print 'Output sidelength (kpc): ' , sideLength
    
    
    armMinPerPixel = float(SideLengthArcMin)/float(kRange)                    
    semiMin = semiMin/armMinPerPixel/2.0
    semiMaj = semiMaj/armMinPerPixel/2.0
     
    
    kpcPerPixelOut = float(sideLength)/float(kRange)         # Output has square aspect ratio
    cmPerPixelOut  = kpcPerPixelOut*kpc2cm                   # Output cm/pixel
    solidAngle     = 1/((dist*kpc2cm)**2.0)              # 4 pi * ratio of 1 cm^2 to the entire sphere given distance
    
    srPerPixel = (kpcPerPixelOut/(dist))**2   # 4 pi * ratio of areas
    SqArcSecPerPixel = sqArcSecPerSr*srPerPixel     # ArcSec^2/Pixel
    pixelPerArcSecSq  = 1/(SqArcSecPerPixel)   # pixel/ArcSec^2  (Ratio of areas times 4 pi * conversion)
    beamConversion    = pixelPerArcSecSq * (math.pi * (hpbw/2)**2) / math.log(2)         # pixel/ArcSec^2  * arcSec^2/beam   
    
    image = numpy.zeros((kRange, kRange))
    
    posAng = -posAng/180*math.pi # inverse rotation and convert to rad.
    flux = 0.0                         # total flux sum
    for i in range(0,kRange):                  # x loop
        for j in range(0,kRange):              # y loop
            
                x = float(i-kRange/2.0)
                y = float(j-kRange/2.0)
                if ( ((x*math.cos(posAng) + y*math.sin(posAng))**2.0/semiMin**2.0 + (x* math.sin(posAng)-y*math.cos(posAng))**2.0/semiMaj**2.0)<=1):
                    image[i,j]= scidata[i,j]  
                    flux += scidata[i, j]         
    
    print 'flux = ', flux/beamConversion, ' mJy, Distance = ' , dist
    
    
#    import matplotlib.pyplot as plt  #@UnresolvedImport
#    import matplotlib.image as mpimg #@UnresolvedImport
#    import matplotlib.cm as cm
#    from matplotlib.patches import Ellipse
#    import matplotlib.patches as mpatches
#    import matplotlib
#    from matplotlib.collections import PatchCollection
#    
#    imgplot = plt.imshow(image)
#    
#    plt.show()
    
    
    
    
#    
#print '\nm31\n'
#intEllipse('m31_out_DM_ONLY.fits',190.0,60.0,700.0,-55,48.)
#
#print '\n2683\n'
#intEllipse('ngc2683_out_DM_ONLY.fits',9.3,2.2,10182.0,-46.5,48.)
#intEllipse('ngc2683_out_DM_ONLY.fits',9.3,2.2,7959.0,-46.5,48.)
#intEllipse('ngc2683_out_DM_ONLY.fits',9.3,2.2,12405.0,-46.5,48.)
#
#print '\n4448\n'
#intEllipse('ngc4448_out_DM_ONLY.fits',3.9,1.4,13000.0,7.9,60.)
#intEllipse('ngc4448_out_DM_ONLY.fits',3.9,1.4,9700.0,7.9,60.)
#intEllipse('ngc4448_out_DM_ONLY.fits',3.9,1.4,47400.0,7.9,60.)
#
#print '\n4866\n'
#intEllipse('ngc4866_out_DM_ONLY.fits',6.3,1.3,21900.0,7.9,54.)
#intEllipse('ngc4866_out_DM_ONLY.fits',6.3,1.3,16000.0,7.9,54.)
#intEllipse('ngc4866_out_DM_ONLY.fits',6.3,1.3,29500.0,7.9,54.)
#
#print '\n1350\n'
#intEllipse('ngc1350_out_DM_ONLY.fits',5.2,2.8,20938.0,90.0,48.)
#
#print '\n7814\n'
#intEllipse('ngc7814_out_DM_ONLY.fits',5.5,2.3,17171.0,45.0,48.)
#
#print '\n4394\n'
#intEllipse('ngc4394_out_DM_ONLY.fits',3.6,3.2,16800.0,50.0,54.)
#print '\n4698\n'
#intEllipse('ngc4698_out_DM_ONLY.fits',4.0,2.5,23650.0,80.0,54.)
#intEllipse('ngc4698_out_DM_ONLY.fits',4.0,2.5,16909.0,80.0,54.)
#intEllipse('ngc4698_out_DM_ONLY.fits',4.0,2.5,30391.0,80.0,54.)

# m31
#  PARAMETERS:  emissFilename, fileout, r_kpc, z_kpc, inclination, objectDist, frequency, hpbw,contour_n,width
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_28.4_16.0_1e+29_0.3_1.0_1.0_22.0', 'm31_out_DM_ONLY.fits',28.4,32.0,[72.2,-55],700.0, 1.485E9, 48.,-2.0,40.0)
##NGC 2683 # used 8000 initially
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc2683_out_DM_ONLY.fits',20.0,32.0,[82.8,-46.5],10182., 1.485E9, 48.,-3.0,8.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc2683_out_DM_ONLY_closer.fits',20.0,32.0,[82.8,-46.5],7959.0, 1.485E9, 48.,-3.0,8.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc2683_out_DM_ONLY_farther.fits',20.0,32.0,[82.8,-46.5],12405.0, 1.485E9, 48.,-3.0,8.0)
##NGC 4448
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_10.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4448_out_DM_ONLY.fits',10.0,32.0,[71.0,7.9],13000.0, 1.485E9, 60.,-3.0,16.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_10.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4448_out_DM_ONLY_closer.fits',10.0,32.0,[71.0,7.9],9700.0, 1.485E9, 60.0,-3.0,16.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_10.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4448_out_DM_ONLY_farther.fits',10.0,32.0,[71.0,7.9],47400.0, 1.485E9, 60.0,-3.0,16.0)
#
#
#
## NGC 4698  # sigma dist = 6.741 (28%)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4698_out_DM_ONLY.fits',20.0,32.0,[73.44,80.0],23650.0, 1.485E9, 54.0,-3.0,16.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4698_out_DM_ONLY_closer.fits',20.0,32.0,[73.44,80.0],16909.0, 1.485E9, 54.0,-3.0,16.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4698_out_DM_ONLY_farther.fits',20.0,32.0,[73.44,80.0],30391.0, 1.485E9, 54.0,-3.0,16.0)
#
#
## NGC 1350 sigma dist = 3.612 (17%)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_28.4_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc1350_out_DM_ONLY.fits',28.4,32.0,[64.79,90.0],20938.0, 1.485E9, 48.0,-3.0,16.0)
#
## NGC 4394  # Only 1 measurement.    
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_15.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4394_out_DM_ONLY.fits',15.0,32.0,[16.55,50.0],16800.0, 1.485E9, 54.0,-3.0,16.0)
#
## NGC 7814 (22% dist variation)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc7814_out_DM_ONLY.fits',20.0,32.0,[90.,45.0],17171.0, 1.485E9, 48.0,-3.0,16.0)



#===============================================================================
# TESTING
#===============================================================================

#for i in (10,20,30,40):
#    vSyncLow     = 1.02026e8        #These are found in the GALDEF File
#    vSyncFactor  = 1.125            #These are found in the GALDEF File
#    frequency = vSyncLow*vSyncFactor**i
#    calc_emiss('./sync.emiss.extragalactic_DM_Only.fits', 'm31_out_DM_ONLY_test2.fits',20.0,16.0,[0,0],7000.0, frequency, 48.0,-2.0,40.0)


def calc_emiss2(emissFilename, fileout, r_kpc, z_kpc, inclination, objectDist, frequency, hpbw,contour_n,width):
    import pyfits, numpy, math #@UnresolvedImport
    import os, sys
    from sys import stdout
    
    pc2cm = 3.08568025e18
    kpc2cm = 3.08568025e21
    sqArcSecPerSr  = 42545949625.0  # Square arcseconds per sr
    
    # Determine which energy bin to use based on the frequency passed in.
    vSyncLow     = 1.02026e8        #These are found in the GALDEF File
    vSyncFactor  = 1.125            #These are found in the GALDEF File
    
    # Determine correct energy bin
    eBin = int(round(math.log(frequency/vSyncLow)/math.log(vSyncFactor))) 
    print 'Selected Frequency:' + str(frequency/1e9) + 'GHz \nEnergy bin: ' + str(eBin)

    
    # generate a square image.
    #outRes = 200         # Number of pixels for square output image 
    #num_zSteps = 100     # Number of z-steps along the integration line
    outRes=15
    num_zSteps=15
    # load emissivity data
    hdulist = pyfits.open(emissFilename, mode='update')
    scidata = hdulist[0].data
    hdulist.info()
    
    # Read FITS header to get dimensions of emissivity map (pixels)
    dimen = scidata.shape
    iRange = int(dimen[0])
    jRange = int(dimen[1])
    kRange = int(dimen[2])
    lRange = int(dimen[3])
    print "Ranges are", iRange, jRange, kRange, lRange

    if (eBin>=jRange): 
        print "Energy out of range.  Exiting..."
        return
    
    #===========================================================================
    # We now will choose a coordinate system of (x',y',z') with the +z' face of a
    # cube parallel to the observer (output image) so that we integrate along z'
    # axis. Thus we have a cylindrical input map rotated in a cube. This cube is
    # chosen to have side length of 1.5*the diameter of the cylinder.  Points in
    # the primed system are mapped to points in the unprimed system via an inverse
    # euler rotation.  Rotations are unitary so the transpose of the standard
    # rotation matrix yields the inverse.
    #===========================================================================
    
    # Calculate the physical side length of the output map based on object distance and angular size of output. 
    sideLength = 50         # in kpc
    print 'Output sidelength (kpc): ' , sideLength
    
    kpcPerDeltaZ   = float(sideLength)/num_zSteps            # Physical distance for z integration interval
    kpcPerPixelInZ = float(z_kpc)/kRange                     # Input Z has roughly 10:1 aspect ratio z:r
    kpcPerPixelInR = float(r_kpc)/lRange                        
    
    kpcPerPixelOut = float(sideLength)/float(outRes)         # Output has square aspect ratio
    cmPerPixelOut  = kpcPerPixelOut*kpc2cm                   # Output cm/pixel
    solidAngle     = 1/((objectDist*kpc2cm)**2.0)              # 4 pi * ratio of 1 cm^2 to the entire sphere given distance
      
    volPerPixelOut = (cmPerPixelOut**2.0)*kpcPerDeltaZ*kpc2cm  # one integration volume in cm^3
    
    # In initial code was missing this 4*pi 
    srPerPixel = 4*math.pi* (kpcPerPixelOut/objectDist)**2   # 4 pi * ratio of areas
    SqArcSecPerPixel = sqArcSecPerSr*srPerPixel              # ArcSec^2/Pixel
    ## Previous
    #pixelPerArcSecSq  = objectDist**2/((kpcPerPixelOut**2))/sqArcSecPerSr   # pixel/ArcSec^2 
    
    pixelPerArcSecSq  = 1/(SqArcSecPerPixel)   # pixel/ArcSec^2  (Ratio of areas times 4 pi * conversion)
    beamConversion    = pixelPerArcSecSq * (math.pi * (hpbw/2)**2)          # pixel/ArcSec^2  * arcSec^2/beam        
    
    print 'ArcSec^2/pixel: ' + str(1/pixelPerArcSecSq)
    print 'Beam Conversion: ' + str(beamConversion)
    
    image = numpy.zeros((outRes, outRes))
    print 'Beginning Integration'
    
    xs=[]
    ys=[]
    zs=[]
    
    for i in range(0,outRes):                  # x loop
        for j in range(0,outRes):              # y loop
            flux = 0.0                         # total flux sum
            for k in range (0,num_zSteps):     # z-integration
                # Shift coordinates so they are centered on galaxy and convert to kpc before transforming to emiss. coordinates.
                outPosition = ((i-outRes/2.0)*kpcPerPixelOut,(j-outRes/2.0)*kpcPerPixelOut,(k-num_zSteps/2.0)*kpcPerDeltaZ)
                inPosition  = rotate_coordinates(outPosition,inclination[0]/180.0*math.pi,0,inclination[1]/180.0*math.pi)       # Rotate to input emissivity file coords.
                
                r = math.sqrt(inPosition[0]**2.0+inPosition[1]**2.0)                          # Calc r in kpc
                z = inPosition[2]                                                             # Calc z in kpc
                if (abs(z) <= z_kpc/2.0 and r <= float(r_kpc)):                               # Check that we are within the bounds of the input image.
                    xPixel = int(r/float(kpcPerPixelInR))                                     # Calc the actual pixel coordinates
                    yPixel = int(z+z_kpc/2.0)/float(kpcPerPixelInZ)                           # Calc the actual pixel coordinates
                    
                    xs.append(outPosition[0])
                    ys.append(outPosition[1])
                    zs.append(outPosition[2])
                    
                    # Ensure we are have a valid coordinate.  The edge flux may be underestimated since we are truncating, but assuming the emiss file is very small at boundaries this won't be an issue.
                    if (xPixel < lRange and yPixel < kRange):    
                        #density = scidata[0,eBin,yPixel, xPixel]                              # Read density from file
                        density = 1                        
                        flux += volPerPixelOut*density                                        # Volume*energy density
                        
                        
            if(flux!=0):            
                #image[j,i] = flux*frequency*solidAngle # Image in erg/s/cm^2
                image[j,i] = flux # Image in erg/s/cm^2
        
        # Monitor Progress in terminal.
        sys.stdout.write('\r' + str(int(float(i)/outRes*100)))
        sys.stdout.flush()
    
    print '\n'
    print 'Total Volume:' + str(numpy.sum(image))
    print 'Total Flux (Jy):' + str(numpy.sum(image)*1.0e23/frequency) # in Jy
    print 'Complete.  FITS image output to ' + fileout
    
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D#@UnresolvedImport
    import matplotlib.pyplot as plt#@UnresolvedImport
    
    fig = plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    ax = Axes3D(fig)
    ax.scatter(xs,ys,zs)
    ax.set_zlim3d([-25, 25])
    ax.set_xlim3d([-25, 25])
    ax.set_ylim3d([-25, 25])
    #plt.show()

#TESTING PURPOSES
#calc_emiss2('./sync.emiss.extragalactic_DM_Only.fits', 'm31_out_DM_ONLY_test2.fits',20.0,16.0,[25.0,0],7000.0, 1.485E9, 48.0,-2.0,40.0)
