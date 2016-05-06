import math
# Determine which energy bin to use based on the frequency passed in.
vSyncLow     = 1.02026e8        #These are found in the GALDEF File
vSyncFactor  = 1.125            #These are found in the GALDEF File

# Determine correct energy bin

#for i in (22e9,90e9):
#    eBin = int(round(math.log(i/vSyncLow)/math.log(vSyncFactor))) 
#    print 'Selected Frequency:' + str(i/1e9) + 'GHz \nEnergy bin: ' + str(eBin)
    
    
    
#bin = 23
#print 'freq', str(vSyncLow*vSyncFactor**bin)


import total_lum
# Returns flux in Jy from a MW at a given dist in kpc  
def flux(dist, emissfile, diff_width):
    kpc2cm = 3.08568025e21
    
    #integrate(emissfile, diff_width, 16.0, 4.85e9)
    
    
    MW238 = total_lum.integrate(emissfile, diff_width, 32.0, 2.38e9)
    MW149 = total_lum.integrate(emissfile, diff_width, 32.0, 1.49e9) # mw lum at 1.49 GHz in erg/s
    MW485 = total_lum.integrate(emissfile, diff_width, 32.0, 4.85e9)
    MW1500 = total_lum.integrate(emissfile, diff_width, 32.0, 15.0e9)
    e2j149   = 1e23/1.485e9  #
    e2j238   = 1e23/2.38e9  #
    e2j485   = 1e23/4.855e9  #
    e2j1500   = 1e23/4.855e9  #
    e2j500    = 1e23/5.0e9  #
    #flux = MW149*e2j/(4* math.pi*(dist*kpc2cm)**2)
    print 'Flux at 1.49GHz, Dist=' , dist, ': ', MW149*e2j149/(4* math.pi*(dist*kpc2cm)**2), ' Jy'
    print 'Flux at 2.38GHz, Dist=' , dist, ': ', MW238*e2j238/(4* math.pi*(dist*kpc2cm)**2), ' Jy'
    print 'Flux at 4.85GHz, Dist=' , dist, ': ', MW485*e2j485/(4* math.pi*(dist*kpc2cm)**2), ' Jy'
    print 'Flux at 5.00GHz, Dist=' , dist, ': ', MW485*e2j500/(4* math.pi*(dist*kpc2cm)**2), ' Jy'
    print 'Flux at 15.0GHz, Dist=' , dist, ': ', MW1500*e2j1500/(4* math.pi*(dist*kpc2cm)**2), ' Jy'







print '\nm31'
flux(700., './diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_28.4_16.0_1e+29_0.3_1.0_1.0_22.0', 28.4)

print '\nngc 2683'
for i in [10200.,7960.,12400.]:
    flux(i, './diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 20.0)

print '\nngc 4448'
for i in [13000.,9700.,47400.]:
    flux(i, './diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_10.0_16.0_1e+29_0.3_1.0_1.0_22.0', 10.0)
    
print '\nngc 4698'
for i in [23700.,16900.,30400.]:
    flux(i, './diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 20.0)

print '\nngc 7814'
for i in [17200,]:
    flux(i, './diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 20.0)

print '\nngc 1350'
for i in [20900.,]:
    flux(i, './diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_28.4_16.0_1e+29_0.3_1.0_1.0_22.0', 28.4)

print '\nngc 4394'
for i in [16800.,]:
    flux(i, './diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_15.0_16.0_1e+29_0.3_1.0_1.0_22.0', 15.0)



# m31
#  PARAMETERS:  emissFilename, fileout, r_kpc, z_kpc, inclination, objectDist, frequency, hpbw,contour_n,width
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_28.4_16.0_1e+29_0.3_1.0_1.0_22.0', 'm31_out_DM_ONLY.fits',28.4,16.0,[72.2,-55],700.0, 1.485E9, 48.,-2.0,40.0)
#NGC 2683 # used 8000 initially
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc2683_out_DM_ONLY.fits',20.0,16.0,[82.8,-46.5],10182., 1.485E9, 48.,-3.0,8.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc2683_out_DM_ONLY_closer.fits',20.0,16.0,[82.8,-46.5],7959.0, 1.485E9, 48.,-3.0,8.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc2683_out_DM_ONLY_farther.fits',20.0,16.0,[82.8,-46.5],12405.0, 1.485E9, 48.,-3.0,8.0)
##NGC 4448
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_10.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4448_out_DM_ONLY.fits',10.0,16.0,[71.0,7.9],13000.0, 1.485E9, 60.,-3.0,16.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_10.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4448_out_DM_ONLY_closer.fits',10.0,16.0,[71.0,7.9],9700.0, 1.485E9, 60.0,-3.0,16.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_10.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4448_out_DM_ONLY_farther.fits',10.0,16.0,[71.0,7.9],47400.0, 1.485E9, 60.0,-3.0,16.0)
#
#
#
## NGC 4698  # sigma dist = 6.741 (28%)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4698_out_DM_ONLY.fits',20.0,16.0,[73.44,80.0],23650.0, 1.485E9, 54.0,-3.0,16.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4698_out_DM_ONLY_closer.fits',20.0,16.0,[73.44,80.0],16909.0, 1.485E9, 54.0,-3.0,16.0)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4698_out_DM_ONLY_farther.fits',20.0,16.0,[73.44,80.0],30391.0, 1.485E9, 54.0,-3.0,16.0)
#
#
## NGC 1350 sigma dist = 3.612 (17%)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_28.4_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc1350_out_DM_ONLY.fits',28.4,16.0,[64.79,90.0],20938.0, 1.485E9, 48.0,-3.0,16.0)
#
## NGC 4394  # Only 1 measurement.    
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_15.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc4394_out_DM_ONLY.fits',15.0,16.0,[16.55,50.0],16800.0, 1.485E9, 54.0,-3.0,16.0)
#
## NGC 7814 (22% dist variation)
#calc_emiss('./diffsynchrotron_emiss_54_grid.60.0_4.0_1.8_20.0_16.0_1e+29_0.3_1.0_1.0_22.0', 'ngc7814_out_DM_ONLY.fits',20.0,16.0,[90.,45.0],17171.0, 1.485E9, 48.0,-3.0,16.0)






#
#for dist in (47400,16000,29500):
#    #print 'Total Flux at ', dist, 'kpc:', flux(dist), ' Jy'
#    flux(dist)
#    print ''
