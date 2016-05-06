#===============================================================================
# Author: Eric Carlson 
# E-mail: erccarls@ucsc.edu
# Description: Integrates galprop emissivity maps for total luminosity
# Last Modified: May 23 2012
#===============================================================================
from scipy import *
import sys, pickle
sys.path.append('./pyfits/lib')
import pyfits
import os
import pickle, numpy, math
import matplotlib.pyplot as plt  #@UnresolvedImport
import matplotlib #@UnresolvedImport
from operator import itemgetter



#===============================================================================
# Integrates the total flux on a relative scale.  Height it FULL WIDTH
#===============================================================================
def integrate(emissFilename, widthKpc, heightKpc, frequency=1.485E9):
    import numpy as np
    
    # Determine which energy bin to use based on the frequency passed in.
    vSyncLow     = 1.02026e8        #These are found in the GALDEF File
    vSyncFactor  = 1.125            #These are found in the GALDEF File
    
    pc2cm = 3.08568025e18
    #frequency=1.485E9
    #eBin= 23
    eBin = int(round(math.log(frequency/vSyncLow)/math.log(vSyncFactor)))
    
    #jyConvert = 1.0e23/frequency
    # load emissivity data
    hdulist = pyfits.open(emissFilename, mode='update')
    scidata = hdulist[0].data
    #hdulist.info()
    # Read FITS header to get dimensions of emissivity map (pixels)
    dimen = scidata.shape
    iRange = int(dimen[0])
    jRange = int(dimen[1])
    kRange = int(dimen[2])
    lRange = int(dimen[3])
    #print "Ranges are", iRange, jRange, kRange, lRange
    
    lumTotal = 0
    #for q in range(0,jRange):# Total Haze Lum
    #for q in range(0,59):    # WMAP haze Lum
    #    frequency = vSyncLow* vSyncFactor**q
    #    flux = np.float64(0.0) # Total flux 
    #    for i in range(0,lRange):                  # r loop
    #        for j in range(0,kRange):              # h loop
    #            flux += 4 * math.pi* scidata[0,q,j, i]*frequency * (math.pi* heightKpc/kRange * ((i+1)**2-i**2)*(widthKpc/lRange)**2 * (1000*pc2cm)**3) # Density times wedge volume 
    #    print "L[", frequency/1e9," GHz]:", flux
    #    lumTotal += flux
    
    
    #frequency = vSyncLow*vSyncFactor**eBin
    lum = np.float64(0.0)
    for i in range(0,lRange):                  # r loop
        for j in range(0,kRange):              # h loop
            lum += 4 * math.pi* scidata[0,eBin,j, i]*frequency * (math.pi* heightKpc/kRange * ((i+1)**2-i**2)*(widthKpc/lRange)**2 * (1000*pc2cm)**3) # Density times wedge volume* 4 pi  
    #print "L[", frequency/1e9," GHz]:", lum
    lumTotal += lum
    return lumTotal
    

def integrate2(emissFilename, widthKpc, heightKpc, frequency=1.485E9):
    import numpy as np
    
    # Determine which energy bin to use based on the frequency passed in.
    vSyncLow     = 1.02026e8        #These are found in the GALDEF File
    vSyncFactor  = 1.125            #These are found in the GALDEF File
    
    pc2cm = 3.08568025e18
    #frequency=1.485E9
    #eBin= 23
    eBin = int(round(math.log(frequency/vSyncLow)/math.log(vSyncFactor)))
    
    #jyConvert = 1.0e23/frequency
    # load emissivity data
    hdulist = pyfits.open(emissFilename, mode='update')
    scidata = hdulist[0].data
    #hdulist.info()
    # Read FITS header to get dimensions of emissivity map (pixels)
    dimen = scidata.shape
    iRange = int(dimen[0])
    jRange = int(dimen[1])
    kRange = int(dimen[2])
    lRange = int(dimen[3])
    #print "Ranges are", iRange, jRange, kRange, lRange
    
    z_eval = 1.*float(kRange)/float(heightKpc) # 1 kpc
    z_flux=0.0
    yidx = round((float(kRange)/2.0) + z_eval) # 2kpc z
    #print yidx, kRange
    #if yidx>kRange:
    #    return(0,0)
    
    r, r_z = [],[]
    
    # Loop through points at z = 1kpc and sum.  This gives LOS integral at 2kpc for edge on.
    for i in range(lRange):
        z_flux += scidata[0,eBin,yidx, i]
        #print i,yidx,scidata[0,eBin,yidx, i]
        r.append(i)
        r_z.append(scidata[0,eBin,yidx, i])
    # Perform LOS integral in R
    r_eval = 5. * float(lRange)/widthKpc  # point to evaluate along r (number in kpc)
    r_flux = 0.0
    
    for x in range(80): # x is the LOS parameter
        r = round(r_eval/math.cos(math.atan(x/2.0/r_eval))) # the 2.0 is just to increase the resolution.
        if r < lRange:
            r_flux += scidata[0,eBin,round(float(kRange)/2.0), r]

    return (r_flux,z_flux)
    



# for full grid of data and total luminosity
def run(rootDir, gridFilename):
    extension = ""
    files = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "diffsynchrotron_emiss_54_grid" in file) and not "gz" in file]
    print "Found", len(files), "files"
    #for i in files: print i
        
    grid = []   
    newcount = 0;
    # The filenames are encoded with the information for every variable in the order:
    # magnetic field strength (mu_g) _ radial exponential falloff (kpc) _ height exponential falloff (kpc) _ simulation width (kpc) _ simulation height (kpc) _ diffusion constant (cm^2s) _ dark matter local density (GeV/cm^3) _ dark matter density index (alpha)_ (stellar radiation energy dens)
    for file in files:
        parsed = file.split("_")[3:]
        #print parsed
        parsed[0] = parsed[0][5:]
        parsed[7] = parsed[7]#[:3]
        # If these were the earlier runs that didn't have radiation  density, then append a 1
        if (len(parsed)==10):
            parsed.append(25.0)
            parsed.append(.33)
        else:
            #print 'New' + str(parsed)
            newcount += 1
        for i in range(len(parsed)):
            parsed[i] = float64(parsed[i])            
        #print " ".join('%02e' % i for i in parsed)

        flux = float64(integrate(rootDir+file,parsed[3],parsed[4]*2.))
        
        print "Total Flux: " , flux
        
        parsed.append(flux)
        grid.append(parsed)
    
    print 'Total New Files: ' + str(newcount)
    outFile = open(gridFilename, "wb" )
    pickle.dump(grid, outFile)
    
    
def run3(rootDir):
    """
    This run is to look at effects of diffusion height on contours.
    """
    extension = ""
    files = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "diffsynchrotron_emiss_54_grid" in file) and not "gz" in file]
    print "Found", len(files), "files"
    #for i in files: print i
        
    grid = []
    int,r,z = [],[],[]
    newcount = 0;
    
    cont_default = integrate2('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    default = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    # The filenames are encoded with the information for every variable in the order:
    # magnetic field strength (mu_g) _ radial exponential falloff (kpc) _ height exponential falloff (kpc) _ simulation width (kpc) _ simulation height (kpc) _ diffusion constant (cm^2s) _ dark matter local density (GeV/cm^3) _ dark matter density index (alpha)_ (stellar radiation energy dens)
    for file in files:
        parsed = file.split("_")[3:]
        #print parsed
        parsed[0] = parsed[0][5:]
        parsed[7] = parsed[7]#[:3]
        # If these were the earlier runs that didn't have radiation  density, then append a 1
        if (len(parsed)==8):
            parsed.append(1.0)
        else:
            #print 'New' + str(parsed)
            newcount += 1
        for i in range(len(parsed)):
            parsed[i] = float64(parsed[i])
        #print " ".join('%02e' % i for i in parsed)

        flux = float64(integrate(rootDir+file,parsed[3],parsed[4]*2.))/default
        cont = integrate2(rootDir+file,parsed[3],parsed[4]*2.)
        flux_r = cont[0]/cont_default[0]
        flux_z = cont[1]/cont_default[1]

        int.append(flux)
        z.append(flux_z)
        r.append(flux_r)
        
        parsed.append(flux)
        grid.append(parsed[4])
     
    print 'hdiff', grid
    print 'r', r
    print 'z',z
    
    plt.scatter(grid,int,label = 'Integrated',marker='+')
    plt.scatter(grid,r, label = "r=5kpc", marker = 'd',c = 'r')
    plt.scatter(grid,z, label = "z=.75kpc",marker = 's')
    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel(r'$L/L_{MW}$')
    plt.xlabel(r'$h_{diff}$')
    plt.ylim(10**-1,10**1)
    plt.xlim(10**0,10**2)
    plt.show()
    


# For random scan histogram
def run2(rootDir, gridFilename,frequency = 1.485e9,DMOnly = False, contours=False):
    extension = ""
    
    # This one for DM_only files
    files = []
    if (DMOnly==True):
        files = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "diffsynchrotron_emiss_54_grid" in file) and not "gz" in file]
    else:
        # This one for CR only files
        files = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "synchrotron_emiss_54_grid" in file and 'nodm' in file) and not "gz" in file ]
    
    print "Found", len(files), "files"
    #for i in files: print i
    
    count = 0    
    grid = []   
    newcount = 0;
    # The filenames are encoded with the information for every variable in the order:
    # magnetic field strength (mu_g) _ radial exponential falloff (kpc) _ height exponential falloff (kpc) _ simulation width (kpc) _ simulation height (kpc) _ diffusion constant (cm^2s) _ dark matter local density (GeV/cm^3) _ dark matter density index (alpha)_ (stellar radiation energy dens)
    for file in files:
        parsed = file.split("_")[3:]
        #print parsed
        parsed[0] = parsed[0][5:]
        parsed[7] = parsed[7]#[:3]
        # If these were the earlier runs that didn't have radiation  density, then append a 1
        if (len(parsed)==8):
            parsed.append(1.0)
        else:
            #print 'New' + str(parsed)
            newcount += 1
        for i in range(1,len(parsed)):
            parsed[i] = float(parsed[i])
                
        #print " ".join('%02e' % i for i in parsed)
        count+=1
        flux = 0 
        if contours == False:
            flux = integrate(rootDir+file,parsed[3],parsed[4]*2.,frequency = frequency)
        else:
            flux = integrate2(rootDir+file,parsed[3],parsed[4]*2.,frequency = frequency)
        
        
        #Append Total Flux
        grid.append(flux)
        print 'Completed' , count, '/' , len(files)
    
    #print 'Total New Files: ' + str(newcount)
    outFile = open(gridFilename, "wb" )
    pickle.dump(grid, outFile)
    
    
# For random scan histogram
def run2_scatter(rootDir, gridFilename,frequency = 1.485e9,DMOnly = False):
    extension = ""
    
    # This one for DM_only files
    files = []
    if (DMOnly==True):
        files = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "diffsynchrotron_emiss_54_grid" in file) and not "gz" in file]
    else:
        # This one for CR only files
        files = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "synchrotron_emiss_54_grid" in file and 'nodm' in file) and not "gz" in file ]
    
    print "Found", len(files), "files"
    #for i in files: print i
    
    count = 0    
    grid = []   
    newcount = 0;
    # The filenames are encoded with the information for every variable in the order:
    # 0magnetic field strength (mu_g) _ 1radial exponential falloff (kpc) _ 2height exponential falloff (kpc) _ 3simulation width (kpc) _ 4simulation height (kpc) _ 5diffusion constant (cm^2s) _ 6dark matter local density (GeV/cm^3) _ 7dark matter density index (alpha)_ 8(stellar radiation energy dens)
    for file in files[:-1]:
        parsed = file.split("_")[3:]
        #print parsed
        parsed[0] = parsed[0][5:]
        parsed[7] = parsed[7]#[:3]
        # If these were the earlier runs that didn't have radiation  density, then append a 1
        if (len(parsed)==8):
            parsed.append(1.0)
        else:
            #print 'New' + str(parsed)
            newcount += 1
        for i in range(1,len(parsed)):
            parsed[i] = float(parsed[i])
                
        #print " ".join('%02e' % i for i in parsed)
        count+=1
        flux = 0 
        if ( .1<=parsed[6] <=.5 ) and (.8<=parsed[7] < 1.2) and (16<=parsed[9]<=28):
            flux = integrate(rootDir+file,parsed[3],parsed[4]*2.,frequency = frequency)
       
        
        
        #Append Total Flux
        grid.append(flux)
        print 'Completed' , count, '/' , len(files)

    #print 'Total New Files: ' + str(newcount)
    outFile = open(gridFilename, "wb" )
    pickle.dump(grid, outFile)
    
    

# Plots the nuissance parameter variations
def plotParam(gridFilename, fileOut):
    
    grid = pickle.load(open(gridFilename, "r" ))
    mwValues = numpy.array([60.0, 4.0, 1.8,20, 16,1.0e29,0.3,1.0,1.0,22.0,25.0,.33,integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)])
    
    gridoutX = []
    gridoutY = []
    for column in range(len(mwValues)-1):
        values  = []
        grid = sorted(grid,key=itemgetter(column))
        for row in grid:
            if (row[column] != mwValues[column]):
                values.append(row)
                #print " ".join('%02e' % i for i in row)
            
            
            if  numpy.all((row[:-1]==mwValues[:-1])):
                values.append(row)
                #print " ".join('%02e' % i for i in row)
               
        seriesX = []
        seriesY = [] 
        for row in values:
            seriesX.append((row[column])/mwValues[column])
            seriesY.append((row[len(row)-1]/mwValues[len(mwValues)-1]))
        gridoutX.append(seriesX)
        gridoutY.append(seriesY)
        
    plt.clf()
    fig1 = plt.figure(1,figsize=(6,6))
    #fig1.suptitle("Total Luminosity",fontsize=14)
    
    ax1 = fig1.add_subplot(3,1,3)
    #ax1.set_ylabel(r'$L/L_{MW}$',fontsize=12)
    ax1.set_xlabel(r'$x/x_{mw}$',fontsize=12)
    ax1.set_autoscaley_on(False)
    ax1.set_autoscalex_on(False)
    plt.xlim([.05,10])
    plt.ylim([.1,10])
    
    plots1 = []
    for i in [0,1,2]:
        t, = ax1.loglog(gridoutX[i],gridoutY[i], ms = 4.0)
        plots1.append(t)
    plots1.append(ax1.loglog((.1,10),(1,1),'--',color="black"))
    ax2 = fig1.add_subplot(3,1,2,sharex=ax1)
    ax2.set_ylabel(r'$L/L_{MW}$',fontsize=12)
    
    plt.setp( ax2.get_xticklabels(), visible=False)
    #ax2.set_xlabel(r'$x/x_{mw}$',fontsize=12)
    ax2.set_autoscaley_on(False)
    #plt.xlim([.5,2])
    plt.ylim([.1,10])

    plots2 = []
    for i in [5,4,3,8,10,11]:
        t, = ax2.loglog(gridoutX[i],gridoutY[i], ms = 4.0)
        plots2.append(t)
    plots2.append(ax2.loglog((.1,10),(1,1),'--',color="black"))
    ax3 = fig1.add_subplot(3,1,1,sharex=ax1)
    plt.setp( ax3.get_xticklabels(), visible=False)
    #ax3.set_ylabel(r'$L/L_{MW}$',fontsize=12)
    #ax3.set_xlabel(r'$x/x_{mw}$',fontsize=12)
    plt.xlim([.1,10])
    
    
    plots3 = []
    for i in [6,7,9]:
        t, = ax3.loglog(gridoutX[i],gridoutY[i],  ms = 4.0)
        plots3.append(t)
    
    plots3.append(ax3.loglog((.1,10),(1,1),'--',color="black"))
    ax3.set_autoscaley_on(False)
    #plt.xlim([.5,2])
    plt.ylim([.1,10])
    
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')    
    
    
    ax1.legend(plots1,["$B_0$","$r_0$","$z_0$"], ncol=1, fancybox=True, shadow=False,prop=fontP,loc = 3, borderaxespad=.5,labelspacing = 0)
    ax2.legend(plots2,["$D_0$", r"$h_{diff}$", r'$R_{diff}$', r'$u_{rad}$',r'$v_{A}$',r'$\gamma_D$'], loc=8, ncol=6,labelspacing = 0, fancybox=True, shadow=False,borderaxespad=.5,prop=fontP,handleheight=.8)
    ax3.legend(plots3,[r"$\rho_0$", r"$\alpha$",r"$r_s$"], ncol=1, fancybox=True, shadow=False,loc = 3,borderaxespad=.5,labelspacing = 0,prop=fontP)
    
    
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(fileOut + '.pdf')
    
    fig1.savefig(pp, format='pdf')
    #fig2.savefig(str(fileOut+"2")+ '.png')
    
    print "Figures saved to ", str(fileOut)+ '.pdf\n',
    pp.close()
    
    fig1.savefig(fileOut + '.png', format='png')
    
    plt.hold(True)
    plt.show()
    
    
def plotHistScatter(gridFilename,CSVFile, gridCRFilename, fileOut,secondPlot = ''):
    import pickle, numpy, math
    import matplotlib.pyplot as plt  #@UnresolvedImport
    import matplotlib #@UnresolvedImport
    import csv
    
    LUM_MW_DM_ONLY = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    LUM_MW_CR_ONLY = integrate('./synchrotron_emiss_54_grid.nodm.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    LUM_MW_ALL = integrate('./synchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    L_MW_DMONLY_30ghz = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=30.E9)
    L_MW_DMONLY_44ghz = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=44.E9)
    
    
    #LUM_MW_DM_ONLY = L_MW_DMONLY_30ghz
    LUM_MW_DM_ONLY44 = L_MW_DMONLY_44ghz
    
    grid = pickle.load(open(gridFilename, "r" ))
    grid2 = []
    
    gridCR = pickle.load(open(gridCRFilename, "r" ))
    grid2CR = []
    
    
    grid44 = pickle.load(open('grid_random_out_DM_ONLY_44ghz.pickle', "r" ))
    grid244 = []
    gridCR44 = pickle.load(open('grid_random_out_CR_ONLY_44ghz.pickle', "r" ))
    grid2CR44 = []
    
    
    suppressed20 = 0.
    suppressed100 = 0.
    # load DM_grid
    
    for i in grid:
        
        if (math.isnan(i)==False ):
            grid2.append(i/LUM_MW_DM_ONLY)  
            print i/LUM_MW_DM_ONLY
            if i/LUM_MW_DM_ONLY<.05:
                suppressed20+=1.
            if i/LUM_MW_DM_ONLY<.01:
                suppressed100+=1.
    print '30 ghz< .01 ' , float(suppressed100)/float(len(grid2))
    print '30 ghz< .05 ' , float(suppressed20)/float(len(grid2))
    
    
    
    # load CR_ONLY_grid
    for i in gridCR:
        if (math.isnan(i)==False):
            grid2CR.append(i/LUM_MW_DM_ONLY)  
            
    sup44_20 =0.
    sup44_100 =0.
    # load DM_grid
    for i in grid44:
        if (math.isnan(i)==False):
            grid244.append(i/LUM_MW_DM_ONLY44)  
            #print i/LUM_MW_DM_ONLY
            if i/LUM_MW_DM_ONLY44<.05: 
                sup44_20 +=1.
            if i/LUM_MW_DM_ONLY44<.01:
                sup44_100 +=1.
                
    print '44 ghz< .01 ' , float(sup44_100)/float(len(grid44))
    print '44 ghz< .05 ' , float(sup44_20)/float(len(grid44))
    
    # load CR_ONLY_grid
    for i in gridCR44:
        if (math.isnan(i)==False):
            grid2CR44.append(i/LUM_MW_DM_ONLY44)            
    
    # Load CSV values
    values= []
    csvread = csv.reader(open(CSVFile, 'rb'))
    for i in csvread:
        values.append(float(i[0]))    
    
    
    plt.scatter(grid2CR, grid2,s=5)
    plt.yscale('log')
    plt.xscale('log')
    
    plt.plot(numpy.logspace(-1,4,100),numpy.logspace(-2,3,100))
    
    plt.xlabel(r'$L_{CR}/L_{DM-MW}$')
    plt.ylabel(r'$L_{DM}/L_{DM-MW}$')
    
    
    
    
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(fileOut + '.pdf')
    pp2 = PdfPages('/home/carlson/Dropbox/totalLum_scatter' + '.pdf')
    plt.savefig(pp, format='pdf')
    plt.savefig(pp2, format='pdf')
    #fig2.savefig(str(fileOut+"2")+ '.png')
    print "Figures saved to ", str(fileOut)+ '.pdf\n',
    pp.close()
    plt.savefig(fileOut + '.png', format='png')

    plt.show()



def plotHist(gridFilename,CSVFile, gridCRFilename, fileOut,secondPlot = ''):
    import pickle, numpy, math
    import matplotlib.pyplot as plt  #@UnresolvedImport
    import matplotlib #@UnresolvedImport
    import csv
    
    LUM_MW_DM_ONLY = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    LUM_MW_CR_ONLY = integrate('./synchrotron_emiss_54_grid.nodm.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    LUM_MW_ALL = integrate('./synchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    L_MW_DMONLY_30ghz = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=30.E9)
    L_MW_DMONLY_44ghz = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=44.E9)
    
    
    #LUM_MW_DM_ONLY = L_MW_DMONLY_30ghz
    LUM_MW_DM_ONLY44 = L_MW_DMONLY_44ghz
    
    grid = pickle.load(open(gridFilename, "r" ))
    grid2 = []
    
    gridCR = pickle.load(open(gridCRFilename, "r" ))
    grid2CR = []
    
    
    grid44 = pickle.load(open('grid_random_out_DM_ONLY_44ghz.pickle', "r" ))
    grid244 = []
    gridCR44 = pickle.load(open('grid_random_out_CR_ONLY_44ghz.pickle', "r" ))
    grid2CR44 = []
    
    
    suppressed20 = 0.
    suppressed100 = 0.
    # load DM_grid
    for i in grid:
        if (math.isnan(i)==False):
            grid2.append(i/LUM_MW_DM_ONLY)  
            #print i/LUM_MW_DM_ONLY
            if i/LUM_MW_DM_ONLY<.05:
                suppressed20+=1.
            if i/LUM_MW_DM_ONLY<.01:
                suppressed100+=1.
    print '30 ghz< .01 ' , float(suppressed100)/float(len(grid2))
    print '30 ghz< .05 ' , float(suppressed20)/float(len(grid2))
    
    
    
    # load CR_ONLY_grid
    for i in gridCR:
        if (math.isnan(i)==False):
            grid2CR.append(i/LUM_MW_DM_ONLY)  
            
    sup44_20 =0.
    sup44_100 =0.
    # load DM_grid
    for i in grid44:
        if (math.isnan(i)==False):
            grid244.append(i/LUM_MW_DM_ONLY44)  
            #print i/LUM_MW_DM_ONLY
            if i/LUM_MW_DM_ONLY44<.05: 
                sup44_20 +=1.
            if i/LUM_MW_DM_ONLY44<.01:
                sup44_100 +=1.
                
    print '44 ghz< .01 ' , float(sup44_100)/float(len(grid44))
    print '44 ghz< .05 ' , float(sup44_20)/float(len(grid44))
    
    # load CR_ONLY_grid
    for i in gridCR44:
        if (math.isnan(i)==False):
            grid2CR44.append(i/LUM_MW_DM_ONLY44)            
    
    # Load CSV values
    values= []
    csvread = csv.reader(open(CSVFile, 'rb'))
    for i in csvread:
        values.append(float(i[0]))    
    
    s = plt.hist(values,bins=10**numpy.linspace(-3, 6,num=45),color='r')
    r = plt.hist(grid2,bins=10**numpy.linspace(-3, 6,num=45),color='b')
    rCR = plt.hist(grid2CR,bins=10**numpy.linspace(-3, 6,num=45),color='g')
    r44 = plt.hist(grid244,bins=10**numpy.linspace(-3, 6,num=45),color='r')
    rCR44 = plt.hist(grid2CR44,bins=10**numpy.linspace(-3, 6,num=45),ls='--', color='r')

    sel = plt.hist(numpy.sort(values)[:7],bins=10**numpy.linspace(-3, 6,num=45),ls='--', color='k')
    bins=10**numpy.linspace(-3, 6,num=45)
    width = [bins[i+1]-bins[i] for i in range(len(bins)-1)]
    
    plt.clf()
    #plt.bar(sel[1][:-1], sel[0]/float(len(values)),color='g',width = width)
    
    #plt.step(s[1][1:], s[0]/float(len(values)),color='g')
    plt.step(r[1][1:], r[0]/float(len(grid2)),ls='-', color='b')
    plt.step(rCR[1][1:], rCR[0]/float(len(grid2CR)),ls=':', color='b')
    plt.step(r44[1][1:], r44[0]/float(len(grid244)),ls='-', color='r')
    plt.step(rCR44[1][1:], rCR44[0]/float(len(grid2CR44)),ls=':', color='r')
    plt.xscale('log')
    #plt.axvline(LUM_MW_ALL/LUM_MW_DM_ONLY,color='k')
    #plt.text(LUM_MW_ALL/LUM_MW_DM_ONLY+3, .063, r'$MW_{DM+CR}$', rotation=90.0)
    plt.ylabel(r'$f$')
    #plt.xlim(10**-3,10**5)
    plt.xlabel(r'$L/L_{MW,DM}$',fontsize=14)
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    
    #plt.legend(['Condon Atlas','DM Only','CR Only'],loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.)
    plt.legend(['30GHz DM','30GHz CR','44GHz DM','44 GHz CR'],loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.)

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(fileOut + '.pdf')
    pp2 = PdfPages('/home/carlson/Dropbox/totalLum' + '.pdf')
    plt.savefig(pp, format='pdf')
    plt.savefig(pp2, format='pdf')
    #fig2.savefig(str(fileOut+"2")+ '.png')
    print "Figures saved to ", str(fileOut)+ '.pdf\n',
    pp.close()
    plt.savefig(fileOut + '.png', format='png')

    plt.show()
    
    
def plotHistContour(gridFilename,CSVFile, gridCRFilename, fileOut,secondPlot = ''):
    import pickle, numpy, math
    import matplotlib.pyplot as plt  #@UnresolvedImport
    import matplotlib #@UnresolvedImport
    import csv
    
    LUM_MW_DM_ONLY = integrate2('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    LUM_MW_CR_ONLY = integrate2('./synchrotron_emiss_54_grid.nodm.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    LUM_MW_ALL = integrate2('./synchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    
    LUM_MW_DM_ONLY_total = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=1.49E9)
    
    grid = pickle.load(open(gridFilename, "r" ))
    grid2 = [[],[]]
    
    gridCR = pickle.load(open(gridCRFilename, "r" ))
    grid2CR = [[],[]]
    
    
    grid_total = pickle.load(open('./grid_random_out.pickle', "r" ))
    grid2_total = []
    for i in grid_total:
        if (math.isnan(i)==False):
            grid2_total.append(i/LUM_MW_DM_ONLY_total)  
    
    suppressedr20 = 0.
    suppressedr100 = 0.
    suppressedz20 = 0.
    suppressedz100 = 0.
    # load DM_grid
    for i in grid:
        if (math.isnan(i[0])==False and math.isnan(i[1])==False):
            grid2[0].append(i[0]/LUM_MW_DM_ONLY[0])
            grid2[1].append(i[1]/LUM_MW_DM_ONLY[1])  
            #print i/LUM_MW_DM_ONLY
            if i[0]/LUM_MW_DM_ONLY[0]<.05:
                suppressedr20+=1.
            if i[0]/LUM_MW_DM_ONLY[0]<.01:
                suppressedr100+=1.
            if i[1]/LUM_MW_DM_ONLY[1]<.05:
                suppressedz20 += 1.
            if i[1]/LUM_MW_DM_ONLY[1]<.01:
                suppressedz100 += 1.
    print '1.49 ghz (r,z)< .01 ' , float(suppressedr100)/float(len(grid2[0])),float(suppressedz100)/float(len(grid2[1]))
    print '1.49 ghz (r,z)< .05 ' , float(suppressedr20)/float(len(grid2[0])),float(suppressedz20)/float(len(grid2[1]))
    
    
    
    # load CR_ONLY_grid
    for i in gridCR:
        if (math.isnan(i[0])==False and math.isnan(i[1])==False):
            grid2CR[0].append(i[0]/LUM_MW_DM_ONLY[0])
            grid2CR[1].append(i[1]/LUM_MW_DM_ONLY[1])    
    
    # Load CSV values
    values= []
        
    
    
    r_total = plt.hist(grid2_total,bins=10**numpy.linspace(-3, 4,num=45),color='b')
    
    r = plt.hist(grid2[0],bins=10**numpy.linspace(-3, 4,num=45),color='b')
    rCR = plt.hist(grid2CR[0],bins=10**numpy.linspace(-3, 4,num=45),color='g')
    
    z = plt.hist(grid2[1],bins=10**numpy.linspace(-3, 4,num=45),color='b')
    zCR = plt.hist(grid2CR[1],bins=10**numpy.linspace(-3, 4,num=45),color='g')
    
    plt.clf()
    
    plt.step(r[1][1:], r[0]/float(len(grid2[0])),ls=':', color='r',label = 'r=5 kpc DM Only')
    #plt.step(rCR[1][1:], rCR[0]/float(len(grid2CR[0])),ls='--', color='r',label = 'r=5 kpc CR Only')
    
    plt.step(z[1][1:], z[0]/float(len(grid2[1])),ls='--', color='b',label = 'z=1 kpc DM Only')
    #plt.step(zCR[1][1:], zCR[0]/float(len(grid2CR[1])),ls='--', color='b',label = 'z=2 kpc CR Only')
    
    plt.step(r_total[1][1:], r_total[0]/float(len(grid2_total)),ls='-', color='g',label = 'Integrated DM Only')
    
    plt.xscale('log')
    plt.ylabel(r'$f$')
    plt.xlabel(r'$L/L_{MW,DM}$',fontsize=14)
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    
    #plt.legend(['r=5 kpc DM Only','r=5 kpc CR Only','z=2 kpc DM Only','z=2 kpc CR Only'],loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.)
    plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.)

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(fileOut + '.pdf')
    plt.savefig(pp, format='pdf')
    #fig2.savefig(str(fileOut+"2")+ '.png')
    print "Figures saved to ", str(fileOut)+ '.pdf\n',
    pp.close()
    plt.savefig(fileOut + '.png', format='png')

    plt.show()
    
    
    
def plotEpsilon(gridFilename, fileOut):    
    
    grid = pickle.load(open(gridFilename, "r" ))
    mwValues = numpy.array([60.0, 4.0, 1.8,20, 16,1.0e29,0.3,1.0,1.0,2.497391e36])
    
    #0 magnetic field strength (mu_g) _ 1radial exponential falloff (kpc) _ 2height exponential falloff (kpc) _ 3simulation width (kpc) _ 4simulation height (kpc) _ 5diffusion constant (cm^2s) _ 6dark matter local density (GeV/cm^3) _ 7dark matter density index (alpha)_ 8(stellar radiation energy dens)
    
    gridoutX = []
    gridoutY = []
    
    Rdiff = mwValues[3]
    hdiff = mwValues[4]
    B0 = mwValues[0]
    z0 = mwValues[2]
    R0 = mwValues[1]
    
    epsBMW =  2* math.pi* B0**2 * z0 * R0**2. *(1-math.exp(-hdiff/z0))*(1-math.exp(-Rdiff/R0)*(1+Rdiff/R0))
    
    for i in grid:
        if (math.isnan(i[9])==False):
            
            Rdiff = i[3]
            hdiff = i[4]
            B0 = i[0]
            z0 = i[2]
            R0 = i[1]
            
            epsB = 2* math.pi* B0**2. * z0 * R0**2 *(1-math.exp(-hdiff/z0))*(1-math.exp(-Rdiff/R0)*(1+Rdiff/R0))
            
            L = i[9]/1.947349e-36
            gridoutX.append(epsB/epsBMW)
            gridoutY.append(L)        
    plt.clf()
    plt.scatter(gridoutX, gridoutY, c ='black', s=1,lw = 1)

    
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$L/L_{MW}$',fontsize=14)
    plt.xlabel(r'$\epsilon_B/\epsilon_{B_{MW}}$',fontsize=14)

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(fileOut + '.pdf')
    plt.savefig(pp, format='pdf')
    plt.savefig(str(fileOut+"2")+ '.png')
    print "Figures saved to ", str(fileOut)+ '.pdf\n',
    pp.close()
    plt.savefig(fileOut + '.png', format='png')

    #plt.show()    


def plotB2(gridFilename, fileOut):    
    
    grid = pickle.load(open(gridFilename, "r" ))
    mwValues = numpy.array([60.0, 4.0, 1.8,20, 16,1.0e29,0.3,1.0,1.0,2.497391e36])
    
    #0 magnetic field strength (mu_g) _ 1radial exponential falloff (kpc) _ 2height exponential falloff (kpc) _ 3simulation width (kpc) _ 4simulation height (kpc) _ 5diffusion constant (cm^2s) _ 6dark matter local density (GeV/cm^3) _ 7dark matter density index (alpha)_ 8(stellar radiation energy dens)
    
    gridoutX = []
    gridoutY = []
    
    Rdiff = mwValues[3]
    hdiff = mwValues[4]
    B0 = mwValues[0]
    z0 = mwValues[2]
    R0 = mwValues[1]
    
    epsBMW =  B0**2
    
    for i in grid:
        if (math.isnan(i[9])==False):
            
            Rdiff = i[3]
            hdiff = i[4]
            B0 = i[0]
            z0 = i[2]
            R0 = i[1]
            
            epsB = B0**2.
            
            L = i[9]/1.947349e-36
            gridoutX.append(epsB/epsBMW)
            gridoutY.append(L)        
    
    plt.clf()
    plt.scatter(gridoutX, gridoutY, c ='black', s=1,lw = 1)

    
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$L/L_{MW}$',fontsize=14)
    plt.xlabel(r'$\epsilon_B/\epsilon_{B_{MW}}$',fontsize=14)

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(fileOut + '.pdf')
    plt.savefig(pp, format='pdf')
    plt.savefig(str(fileOut+"2")+ '.png')
    print "Figures saved to ", str(fileOut)+ '.pdf\n',
    pp.close()
    plt.savefig(fileOut + '.png', format='png')

    #plt.show()    




def plotMass(gridFilename, fileOut):
    
    grid = pickle.load(open(gridFilename, "r" ))
    mwValues = numpy.array([60.0, 4.0, 1.8,20, 16,1.0e29,0.3,1.0,1.0,2.497391e36])
    
    #0 magnetic field strength (mu_g) _ 1radial exponential falloff (kpc) _ 2height exponential falloff (kpc) _ 3simulation width (kpc) _ 4simulation height (kpc) _ 5diffusion constant (cm^2s) _ 6dark matter local density (GeV/cm^3) _ 7dark matter density index (alpha)_ 8(stellar radiation energy dens)
    
    gridoutX1 = []
    gridoutY1 = []
    gridoutX2 = []
    gridoutY2 = []
    gridoutX3 = []
    gridoutY3 = []
    
    
    plt.clf()
    #massMW =  math.pow(8500,mwValues[7])*mwValues[6]*math.log(r_cutoff) # integrating from 1 to r_cutoff
    massMW = mwValues[6]**2
    for i in grid:
        if (math.isnan(i[9])==False):
            alpha = i[7]             # 
            rho = i[6] # relative density compared to MW
            mass = 0
            r0 = 1
            
            L = i[9]/2.497391e36
            if (alpha < 1):
                #mass = rho/(1-alpha)*math.pow(8500,alpha)*(r_cutoff**(1-alpha)-r0**(1-alpha))  / massMW
                mass = rho**2 /(3-2 *alpha)/massMW
                gridoutX1.append(mass)
                gridoutY1.append(L)        
            elif (alpha > 1):
                #mass = rho/(1-alpha)*math.pow(8500,alpha)*(r_cutoff**(1-alpha)-r0**(1-alpha))  / massMW
                mass = rho**2 /(3-2 *alpha)/massMW
                gridoutX2.append(mass)
                gridoutY2.append(L)
    plt.clf()
    p1 = plt.scatter(gridoutX1, gridoutY1, c ='red', s=4,lw = 0,marker = 'o')
    p2 = plt.scatter(gridoutX2, gridoutY2, c ='black', s=4,lw = 0,marker = 'o')
    #plt.autoscale(False, 'x', None)
    #plt.xlim(10**-2,10**3)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$L/L_{MW}$',fontsize=14)
    plt.xlabel(r'$R/R_{MW}$',fontsize=14)
        
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small') 
    plt.legend((p1,p2),(r'$\alpha<1$',r'$\alpha>1$'),loc=2, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.)
    
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(fileOut + '.pdf')
    plt.savefig(pp, format='pdf')
    pp.close()
    
    plt.savefig(fileOut + '.png', format='png')
    
    print "Figures saved to ", str(fileOut)+ '.pdf\n',

    #plt.show()    




def plotDMMass():
    
    vSyncLow     = 1.02026e8        #These are found in the GALDEF File
    vSyncFactor  = 1.125            #These are found in the GALDEF File
    
    rootDir = "/home/carlson/data/mass/"
    extension = ""
    files = [file for file in os.listdir(rootDir) if (file.lower().endswith(extension) and "synchrotron_emiss_54_grid" in file) and not "nodm" in file]
    print "Found", len(files), "files"
    
    flux, mass = [], []    
    
    LUM_MW = integrate('./synchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=22.e9)
    LUM_MW_DM_ONLY = integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=22.e9)
    print LUM_MW_DM_ONLY
    
    X = numpy.logspace(7.5, 11,50)
    norm = LUM_MW_DM_ONLY/(22e9**-1.55)
    normCR = LUM_MW/(22e9**-1.55)
    
    Y = [norm*(i**-2.55) for i in X]
    Y2 = [normCR*(i**-2.55) for i in X]
    
    #plt.plot(X,Y,ls = '-.',c='k')
    plt.plot(X,Y2,label=r'$\beta_H=-2.55$',ls = '-.',c = 'k')
    
    
    masses = [0,3,6,10,14]
    c = ['r','b','g','m','y']
    
    flux149 = []
    
    for i in range(len(masses)):
        mass = masses[i]
        # Determine normalization
        norm = LUM_MW_DM_ONLY/float64(integrate(rootDir+"diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0_25.0_0.33_" + str(mass),20.,32.,frequency=22.e9))
        normCR = LUM_MW/float64(integrate(rootDir+"synchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0_25.0_0.33_" + str(mass),20.,32.,frequency=22.e9))
        freq,flux = [],[]
        fluxCR = []
        
        flux149.append(norm*float64(integrate(rootDir+"diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0_25.0_0.33_" + str(mass),20.,32.,frequency=1.49e9)))
        
        for bin in range(58):
            nu = vSyncLow*vSyncFactor**bin
            freq.append(nu)
            lum = float64(integrate(rootDir+"diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0_25.0_0.33_" + str(mass),20.,32.,frequency=nu))
            flux.append(norm*lum/nu)
        
            lum = float64(integrate(rootDir+"synchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0_25.0_0.33_" + str(mass),20.,32.,frequency=nu))
            fluxCR.append(normCR*lum/nu)
        
        plt.plot(freq,flux,label = r'$m_\chi = $' + str(int(10*10**(.2*mass))) + ' GeV, DM',ls = '--',c = c[i])
        plt.plot(freq,fluxCR,label = r'$m_\chi = $' + str(int(10*10**(.2*mass))) + ' GeV CR+DM',c = c[i])
        
        #print numpy.array(freq) 
        #break
    
    plt.axvline(1.49e9, 0, 1e50,c='k')
    plt.axvline(22e9, 0, 1e50,c='k')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r'$\nu$ [Hz]')
    plt.ylabel(r'L($\nu$)/$\nu$ [erg s$^{-1}$ Hz$^{-1}$]')
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small') 
    
    plt.legend(ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.)
    
    fig = plt.figure(2)
    def f149(mass):
        norm = LUM_MW_DM_ONLY/float64(integrate(rootDir+"diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0_25.0_0.33_" + str(mass),20.,32.,frequency=22.e9))
        
        return norm* integrate(rootDir+"diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0_25.0_0.33_" + str(mass),20.,32.,frequency = 1.49e9)
    
    flux149 = [f149(mass)/1.49e9 for mass in range(14)]
    masses = [ 10*10**(.2*mass) for mass in range(14)]
    
    print masses, flux149
    plt.plot(masses,flux149)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r'$m_\chi$ [GeV]')
    plt.ylabel(r'L($\nu$)/$\nu$ [erg s$^{-1}$ Hz$^{-1}$]')
    plt.show()
    
    
#run3("/home/carlson/diff_height/")
   
# Parameters
# This runs the file and outputs the list for input to the plot routine.    
#run("/sdata/galprop_eric/grid_one_dim/", "grid_out.pickle")
#run("/home/carlson/data/", "grid_out.pickle")
#plotParam("grid_out.pickle","totalLum")

# Histogram
#run2("/sdata/galprop_eric/grid_new/", "grid_random_out.pickle",DMOnly=True)
# CR_ONLY_data 
#run2("/sdata/galprop_eric/grid_new/", "grid_random_out_CR_ONLY.pickle")


#Scatter Correlation Histogram
#run2_scatter("/sdata/galprop_eric/grid_new/", "grid_random_out_scatter.pickle",DMOnly=True)
# CR_ONLY_data 
#run2_scatter("/sdata/galprop_eric/grid_new/", "grid_random_out_CR_ONLY_scatter.pickle")


#run2("/sdata/galprop_eric/grid_new/", "grid_random_out_CR_ONLY_30ghz.pickle",frequency = 3.0e10)
#run2("/sdata/galprop_eric/grid_new/", "grid_random_out_CR_ONLY_44ghz.pickle",frequency = 4.4e10)
#run2("/sdata/galprop_eric/grid_new/", "grid_random_out_DM_ONLY_30ghz.pickle",frequency = 3.0e10,DMOnly=True)
#run2("/sdata/galprop_eric/grid_new/", "grid_random_out_DM_ONLY_44ghz.pickle",frequency = 4.4e10,DMOnly=True)


#plotHist("grid_random_out_DM_ONLY_30ghz.pickle",'hist_obs.csv',"grid_random_out_CR_ONLY_30ghz.pickle","scanHist_30ghz")
#plotHist("grid_random_out.pickle",'hist_obs.csv',"grid_random_out_CR_ONLY.pickle","scanHist")
plotHistScatter("grid_random_out.pickle",'hist_obs.csv',"grid_random_out_CR_ONLY.pickle","scanHist")

# Contour Plots
#run2("/sdata/galprop_eric/grid_new/", "grid_random_out_contour.pickle",DMOnly=True,contours = True)
#run2("/sdata/galprop_eric/grid_new/", "grid_random_out_CR_ONLY_contour.pickle",contours = True)
#plotHistContour("grid_random_out_contour.pickle",'hist_obs.csv',"grid_random_out_CR_ONLY_contour.pickle","scanHist_contour")


#run("/sdata/galprop_eric/haze_random/", "grid_random_out_full.pickle")
#plotEpsilon('grid_random_out_full.pickle', 'epsilonCorrelation')
#plotB2('grid_random_out_full.pickle', 'BCorrelation')
#plotMass('grid_random_out_full.pickle', 'rateCorrelation')

#print 'Total Luminosity: ' + str(integrate('./sync.emiss.extragalactic_DM_Only.fits',20.,32.,frequency=23.E9)) + ' erg/s'
#print 'Total Luminosity: ' + str(integrate('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0',20.,32.,frequency=23.E9)) + ' erg/s'



# CR_ONLY
#print 'Total Luminosity: ' + str(integrate('./synchrotron_emiss_54_grid.nodm.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0',20.,32.)) + ' erg/s'
# ALL
#print 'Total Luminosity: ' + str(integrate('./synchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0',20.,32.,frequency=30.E9)) + ' erg/s'
#print 'Total Luminosity: ' + str(integrate('./synchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0',20.,32.,frequency=44.E9)) + ' erg/s'
#print 'Total Luminosity: ' + str(integrate('./sdata/galprop_eric/haze_random/diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0',20.,32.)) + ' erg/s'


#plotDMMass()



