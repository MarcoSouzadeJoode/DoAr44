#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
import pyshellspec3
import pandas as pd

def read_lc(lfile):
    """
    """
    data = np.loadtxt(lfile, dtype=str)

    names = data[:,0]
    pbands = data[:,1]
    errs = data[:,2].astype('float64')

    return names, pbands, errs

def read_if(ifile):
    """

    :param ifile:
    :return:
    """
    return np.loadtxt(ifile, dtype=str)


ra = (16 + 31 / 60. + 33.4634997336 / 3600.)/24.*360.
dec = -(24 + 27 / 60. + 37.158728076 / 3600.)
dir_data = '../data_synthetic_alma'
print(ra, dec)



def main():
    """
    Main loop
    :return:
    """
    # read the data
    obs = []

   
    """
    # ALMA data
    IF_data_filename = '../data/VIS_cont_ABC.dat'
    obs.append(pyshellspec3.IFData(filename=IF_data_filename, location='alma', ra=ra, dec=dec, format='ascii'))
    """
    
    # VLTI data
    IF_data_filename = '../data/VLTI_multiplied.dat'
    obs.append(pyshellspec3.IFData(filename=IF_data_filename, location='VLT', ra=ra, dec=dec, format='ascii'))
    
    #HALFA
    HALFA_combined = '../data/HALFA_combined_NEW'
    obs.append(pyshellspec3.SPEData(filename=HALFA_combined))
    
    
    #SED data
    SED_data_filename = '../data/SED_NEBULA_DR_cgs_flux_inner'
    obs.append(pyshellspec3.SEDData(filename=SED_data_filename))
    
    # VLTI data
    CLO_data_filename = '../data/Clo_DoAr44'
    obs.append(pyshellspec3.IFData(filename=CLO_data_filename, location='VLT', ra=ra, dec=dec, format='ascii'))
    
    
    #SED data
    #SED_data_filename = '../data_synthetic_alma/SED_NEBULA_DR'
    #obs.append(pyshellspec3.SEDData(filename=SED_data_filename))
    
    
    #UVES_halpha
    #UVES_HALPHA_filename = "../data_synthetic_alma/UVES_HALFA"
    #obs.append(pyshellspec3.SPEData(filename=UVES_HALPHA_filename))
    
    # XSHOOTER_HALFA
    #XSHOOTER_HALFA = '../data_synthetic_alma/XSHOOTER_HALFA'
    #obs.append(pyshellspec3.SPEData(filename=XSHOOTER_HALFA))
    
    # one line
    #halpha = '../data_synthetic_alma/just_halpha'
    #obs.append(pyshellspec3.SPEData(filename=halpha))

    
    # XSHOOTER_PASCHEN
    #XSHOOTER_PASCHEN = '../data_synthetic_alma/XSHOOTER_PASCHEN'
    #obs.append(pyshellspec3.SPEData(filename=XSHOOTER_PASCHEN))
    """
    
    # UVES_HBETA
    UVES_HBETA = '../data_synthetic_alma/UVES_HBETA'
    obs.append(pyshellspec3.SPEData(filename=UVES_HBETA))
    
    # XSHOOTER_PASCHEN_BETA
    XSHOOTER_PASCHEN_BETA = '../data_synthetic_alma/XSHOOTER_PASCHEN_BETA'
    obs.append(pyshellspec3.SPEData(filename=XSHOOTER_PASCHEN_BETA))
    
    # XSHOOTER_PASCHEN_GAMMA, also 1083 helium line
    XSHOOTER_PASCHEN_GAMMA= '../data_synthetic_alma/XSHOOTER_PASCHEN_GAMMA'
    obs.append(pyshellspec3.SPEData(filename=XSHOOTER_PASCHEN_GAMMA))
    """
    

    # construct data class
    data = pyshellspec3.Data(obs)

    # construct the model
    central = pyshellspec3.CentralObject()
    companion = pyshellspec3.Companion()
    nebula = pyshellspec3.Nebula()
    disk = pyshellspec3.Disk()
    jet = pyshellspec3.Jet()
    spot = pyshellspec3.Spot()
    # envelope = pyshellspec3.Envelope()
    shell = pyshellspec3.Shell()
    orbit = pyshellspec3.Orbit()
    
    objs = [central, companion, nebula,shell,disk, jet,spot,orbit]
    model = pyshellspec3.Model(objects=objs)

    
    
    # construct the Interface
    
    pixels = 60
    print("OK")
    itf = pyshellspec3.Interface(model=model, data=data, ncpu=8,
        image_size=pixels*2+1,
        if_phase_precision=3,
        df_phase_precision=3,
        lc_phase_precision=2,
        sed_phase_precision=2,
        spe_phase_precision=2, 
        if_ew_precision=7,
        df_ew_precision=10,
        lc_ew_precision=9,
        sed_ew_precision=10,
        spe_ew_precision=10,
        shellspec_template="template.in",
        shellspec_abundance="/home/marco/Desktop/ALMA/oct2023/working/fitting/abundances",
        use_offset=True,
        use_differential=True,
        exclude_visphi=False, dry_run=False)

    
    N = 35
    #N = 4
    
    z_aspect = 1.0
    # body-frozen grid, line-of-sight grid
    tmp=N + 0.01
    itf.set_parameter('rmdfx1', value=-tmp)
    itf.set_parameter('rmdfx2', value=+tmp)
    itf.set_parameter('rmdfy1', value=-tmp)
    itf.set_parameter('rmdfy2', value=+tmp)
    itf.set_parameter('rmdfz1', value=-tmp/z_aspect)
    itf.set_parameter('rmdfz2', value=+tmp/z_aspect)

    tmp=N
    itf.set_parameter('rmdx1',  value=-tmp)
    itf.set_parameter('rmdx2',  value=+tmp)
    itf.set_parameter('rmdy1',  value=-tmp)
    itf.set_parameter('rmdy2',  value=+tmp)
    itf.set_parameter('rmdz1',  value=-tmp/z_aspect)
    itf.set_parameter('rmdz2',  value=+tmp/z_aspect)

    tmp=N / pixels
    aspect_sight = z_aspect
    itf.set_parameter('stepfx', value=tmp)
    itf.set_parameter('stepfy', value=tmp)
    itf.set_parameter('stepfz', value=tmp/aspect_sight)

    itf.set_parameter('stepx', value=tmp)
    itf.set_parameter('stepy', value=tmp)
    itf.set_parameter('stepz', value=tmp/aspect_sight)

    # set parameters


    # STAR ---------
    itf.set_parameter('istar', value=1)
    itf.set_parameter('inebl', value=0)  # nebula

    # R. Ricci et al. 2010.
    itf.set_parameter('rstar', value=2)
    
    # modified
    itf.set_parameter('Tstar', value=4750)
    
    # END STAR ----------


    
    
    itf.set_parameter('ichemc', value=1)
    itf.set_parameter('ielnd', value=1)
    itf.set_parameter('ithom', value=1) # scattering
    itf.set_parameter('irayl', value=1)
    
    itf.set_parameter('imie', value=0)
    itf.set_parameter('imiepf', value=0)
    itf.set_parameter('iline', value=1)

    itf.set_parameter('ispot', value=0)
    itf.set_parameter('ijet', value=0)

    
    


    itf.set_parameter('dd'      , value=146.3                , fitted=False, vmin=305., vmax=330.)
    
    # TODO: najít si gamma
    itf.set_parameter('vgamma'  , value=100               , fitted=False, vmin=-100., vmax=100.)
    
    

    # skoro jako testovací částice v centrálním poli
    # 1 rok a 1 au => 1 M_sun, snad!!


    # změřeno v CARTA z Band 6, 7 pozorování. Zde ale nehraje roli
    dinc = 40
    itf.set_parameter('dinc'    , value=dinc                 , fitted=False, vmin=10., vmax=30.)
    
    # nastaveno tak, aby odpovídalo M = 1.4 M_sun. takto koresponduje syntetickým spektrům
    itf.set_parameter('q'       , value=0.01               , fitted=False, vmin=0.0, vmax=0.30)
    itf.set_parameter('period'  , value=365 * 0.84           , fitted=False, vmin=0., vmax=1000.)
    itf.set_parameter('asini'   , value=84.0/np.sin(23.0*np.pi/180)*np.sin(dinc*np.pi/180), fitted=False, vmin=53., vmax=100.)
    
    itf.set_parameter('omega_an'   , value=147, fitted=False, vmin=53., vmax=100.)

    # companion --- parasitic star
    # emission is turned off
    # M Cieza et al 2021 : mass of 1.4 Msun
    
    
    
    itf.set_parameter('icomp', value=0)
    itf.set_parameter('ecc'     , value=0.0                  , fitted=False, vmin=0., vmax=0.5)
    itf.set_parameter('rcp', value=0.1, fitted=False, vmin=0.0, vmax=50)
    itf.set_parameter('tempcp', value=6000)
    # ----------------------------------------

    ionu = 42 # 0.87 mm
    ionu =14 # Halpha cca
    ionu = 1
    ior = 40
    iot = 40
    print("IONU = ", ionu)
    itf.set_parameter('ionu'    , value=ionu  , fitted=False, vmin=0, vmax=10000)
    itf.set_parameter('ior'    , value=ior     , fitted=False, vmin=0, vmax=10000)
    itf.set_parameter('iot'    , value=iot     , fitted=False, vmin=0, vmax=10000)
    # NEBULA -----------------------------------------------
    
    
    
    


    
    

    itf.set_parameter('ielnd', value=1)
    itf.set_parameter('ithom', value=1) # scattering
    itf.set_parameter('irayl', value=1)
    
    # matching the calculated mass of the star !!
    # beware!! this value is now hard coded
    # --> OK, done.



    itf.set_parameter('ielnd', value=1)
    itf.set_parameter('ithom', value=1) # scattering
    itf.set_parameter('irayl', value=1)
    itf.set_parameter('imie', value=0) # dust 
    itf.set_parameter('imiepf', value=0)
    itf.set_parameter('iline', value=1)
    itf.set_parameter('iinvnb', value=1)  # inversion of T
    #itf.set_parameter('ivelnb', value=1)  # radial wind
    itf.set_parameter('ishdnb', value=1)  # shadow
    
    itf.set_parameter('rinnb'   , value=1.9 , fitted=False, vmin=0.0, vmax=20000)
    itf.set_parameter('routnb'  , value=18.0 , fitted=False, vmin=6.0, vmax=30000)
    
    
    
    DENSITY = 1e-8
    METALLICITY = 0.02
    
    # DUST
    # temperature given by BB temperature at 8000 rsun = 37 au
    TEMPERATURE = 6000

    
    itf.set_parameter('idisc', value=0)
    itf.set_parameter('ishell', value=3)
    itf.set_parameter('ijet', value=0)
    itf.set_parameter('ispot', value=0)
    itf.set_parameter('inebl', value=1)

    # NEBULA
    itf.set_parameter('rinnb'   , value=16.5, fitted=False, vmin=0.0, vmax=50)
    itf.set_parameter('routnb'   , value=35, fitted=False, vmin=0.0, vmax=90)
    itf.set_parameter('itnb'   , value=3, fitted=False, vmin=0, vmax=10)
    
    itf.set_parameter('tempnb'  , value=2500.0, fitted=False, vmin=1500., vmax=5000.)
   
    itf.set_parameter('densnb'  , value=1e-8, fitted=False, vmin=1e-15, vmax=5e1)
    itf.set_parameter('edennb'  , value=-0.5, fitted=False, vmin=-3.0, vmax=2)
    itf.set_parameter('etmpnb'  , value=-0.7, fitted=False, vmin=-3.0, vmax=3.0)
    itf.set_parameter('aneb'    , value=1.5, fitted=False, vmin=0.0, vmax=12.0)
    # END NEBULA -------------------------------------------------
 
    
    #SHELL ------
    itf.set_parameter('rinsh',  value=3, vmin=1.8, vmax=10)
    itf.set_parameter('rcsh',  value=2, vmin=1, vmax=10)

    itf.set_parameter('vtrbsh', value=1, vmin = 0, vmax = 100)
    
    itf.set_parameter('routsh', value=9, fitted=False, vmin=2, vmax = 45)
    itf.set_parameter('vsh',    value=330, fitted=True, vmin = 200, vmax = 1000)
    itf.set_parameter('denssh', value=2.8e-13, fitted=True, vmin = 1e-14, vmax = 7e-13)
    itf.set_parameter('tempsh', value=10_000,fitted=True, vmin = 9000, vmax = 30000)
    itf.set_parameter('evelsh', value=1.3,fitted=True, vmin = 1, vmax = 4)
    
    
    #END SHELL ----
    
    

    # compute one/fit
    chi2 = itf.compute_chi2(verbose=False)

    #itf.run_fit(fitter='sp_diff_evol', tol=1e-2, maxiter=100)

    # fitter
    #itf.run_fit(fitter='nlopt_nelder_mead', ftol=1e-6, maxiter=200)


    itf.write_iterations()
    itf.set_model_to_shellspec()
    itf.write_template('final.in')
    itf.write_model()

    print("Note: fit.py ended successfully.")
    #sys.exit(0) 

    files = [] 
    for o in obs:
        files.append(o.get_filename().split('/')[-1])
    for f in list(set(files)):
        itf.plot_comparison(filename=f)



if __name__ == '__main__':
    main()


