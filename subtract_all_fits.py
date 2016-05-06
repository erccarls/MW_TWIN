'''
Created on Apr 22, 2012

@author: carlson
'''
import os, sys
sys.path.append('./pyfits/lib')

def subFits(file1,file2,fileout):
    import pyfits, numpy, math #@UnresolvedImport
    
    hdulist1 = pyfits.open(file1, mode='update')
    scidata1 = hdulist1[0].data
    #hdulist1.info()
    
    hdulist2 = pyfits.open(file2, mode='update')
    scidata2 = hdulist2[0].data
    #hdulist2.info()
        
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


#rootDir = "/sdata/galprop_eric/grid_new/"
#outDir = "/sdata/galprop_eric/grid_new/"
#rootDir = "/home/carlson/data/default_final/"
#outDir = "/home/carlson/data/default_final/"
#rootDir = "/sdata/galprop_eric/grid_one_dim/"
#outDir =  "/sdata/galprop_eric/grid_one_dim/"


rootDir = "/home/carlson/data/mass/"
outDir =  "/home/carlson/data/mass/"


#rootDir = "/home/carlson/diff_height/"
#outDir =  "/home/carlson/diff_height/"

extension = ""
nodm = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "synchrotron_emiss_54_grid.nodm" in file) and not "gz" in file]
full = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "synchrotron_emiss_54_grid" in file) and not "gz" in file and not "nodm" in file]
print "Found", len(nodm), "files"
for i in range(0,len(nodm)):
    print nodm[i][31:]
    print full[i][26:]
    for j in range(len(full)):
        if (nodm[i][31:] == full[j][26:]):
            print "Pair Found:"
            print nodm[i][31:]
            print full[j][26:]
            subFits(rootDir+full[j],rootDir+nodm[i],outDir+"diff"+full[j])
    print "Progress: " + str(i) + "/" + str(len(nodm)) 
    
    


