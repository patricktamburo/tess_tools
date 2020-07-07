import lightkurve as lk
import pickle
import pdb
import matplotlib.pyplot as plt
import numpy as np
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
import tess_stars2px



def coordinate_propagator(ra, dec, parallax, pmra, pmdec, ref_epoch, target_epoch, tic_id):
    '''Authors:
		Patrick Tamburo, Boston University, July 2020
	Purpose:
        Given a target's position, parallax, proper motion, and reference epoch, propagate space motion to target epoch. 
        See https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html for more examples. 
	Inputs:
        ra (float): The target's DECIMAL RA at your reference epoch
        dec (float): The target's DECIMAL Dec at your reference epoch 
        parallax (float): The target's parallax in milliarcsec
        pmra (flaot): The target's RA proper motion in milliarcsec/year
        pmdec (float): The target's Dec proper motion in millarcsec/year
        ref_epoch (float): Decimal year representing the epoch when coordinates, parallax, and proper motions were measured (e.g., 2015.5 for Gaia measurements)
        target_epoch (str): A string representing the epoch you want to propagate motion to (e.g., '2019-09-25')
        tic_id (str): The target's TIC ID (e.g., 'TIC 326109505')
    Outputs:
        new_ra, new_dec (float): The propagated position of the target at target_epoch.
	TODO:
    '''
    c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=Distance(parallax=parallax*u.mas), pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, obstime=Time(ref_epoch, format='decimalyear'))
    target_obstime = Time(target_epoch)
    c_target_epoch = c.apply_space_motion(target_obstime) #Apply space motion to TESS epoch
    new_ra = c_target_epoch.ra.value
    new_dec = c_target_epoch.dec.value
    return new_ra, new_dec

def tess_pixels(tic_id, tess_ra, tess_dec):
    '''Authors:
		Patrick Tamburo, Boston University, July 2020
	Purpose:
        Given a target's position when it was observed with TESS, returns the column and row of the pixel that the target fell on. 
	Inputs:
        tic_id (str): The target's TIC ID (e.g., 'TIC 326109505')
        tess_ra (float): The target's DECIMAL RA at the TESS epoch when this target was observed
        tess_dec (float): The target's DECIMAL Dec at the TESS epoch when this target was observed
    Outputs:
        out_col_pix, out_row_pix (float): The propagated position of the target at target_epoch.
        out_sec, out_cam, out_ccd: sector, camera, and ccd, respectively
	TODO:
    '''
    tid = int(tic_id.split(' ')[-1])
    out_id, out_eclip_long, out_eclip_lat, out_sec, out_cam, out_ccd, out_col_pix, out_row_pix, scinfo = tess_stars2px.tess_stars2px_function_entry(tid, tess_ra, tess_dec)
    return out_col_pix, out_row_pix, out_sec, out_cam, out_ccd
    
def tpf_grabber(tic_id, sectors):
    '''Authors:
		Patrick Tamburo, Boston University, July 2020
	Purpose:
        Given a target's position when it was observed with TESS, returns the column and row of the pixel that the target fell on. 
	Inputs:
        tic_id (str): The target's TIC ID (e.g., 'TIC 326109505')
        sector (array of ints): Sectors during which the target was observed
    Outputs:
        tpfs (list): List of target pixel files. 
	TODO:
    '''

    tpfs = []
    for i in range(len(sectors)):
        tpf = lk.search_targetpixelfile(tic_id, sector=sectors[i]).download()
        if tpf is not None:
            tpfs.append(tpf)
        else:
            #Sometimes, tpf files don't exist for certain sectors, even though the target may have been observed in them...don't know why. 
            #But if that is the case, just skip that sector. 
            print('ERROR: no data found for {} in sector {:2d}, skipping.'.format(tic_id, sectors[i]))
            continue
    return tpfs
    # tpf = lk.search_targetpixelfile(tic_id,sector=16).download()

    # tpf.plot()

    # plt.figure()
    # plt.imshow(tpf.flux[0],origin='lower')
    # target_x = 1608.116
    # target_y = 422.226
    # plt.plot(target_x-1603.5,target_y-418.5,'rx')
    # target_x_pix = int(np.floor(target_x-1603.5+0.5))
    # target_y_pix = int(np.floor(target_y-418.5+0.5))
    # print('Target on pixel ',target_x_pix,target_y_pix)

    # plt.figure(figsize=(10,5))
    # bg_sub_flux = tpf.flux[:,target_y_pix,target_x_pix]-tpf.flux[:,3,3] #Transiting object is in (4,5) in sector 16
    # lc = lk.LightCurve(time=tpf.time,flux=bg_sub_flux)
    # plt.plot(lc.time,lc.flux,'k.')
    # plt.plot(lc.bin(binsize=20).time,lc.bin(binsize=20).flux,'r.')

if __name__ == "__main__":
    ra = 331.44281956953574
    dec = 53.89447777571114
    parallax = 16.710742916905904
    pmra = -447.77974032397765
    pmdec = -298.5485396420709
    ref_epoch = 2015.5
    target_epoch = '2019-09-25' #Approximately the middle of Sector 16 
    tic_id = 'TIC 326109505'

    #Get coordinates in TESS epoch. 
    tess_epoch_ra, tess_epoch_dec = coordinate_propagator(ra, dec, parallax, pmra, pmdec, ref_epoch, target_epoch, tic_id)
    print('')
    print('{:>30s} {}, {}'.format('Original epoch RA, Dec:', np.round(ra,5), np.round(dec,5)))
    print('{:>30s} {}, {}'.format('TESS epoch RA, Dec:', np.round(tess_epoch_ra,5), np.round(tess_epoch_dec,5)))
    print('')

    #Get the pixels it fell on. 
    target_x, target_y, sectors, cameras, ccds = tess_pixels(tic_id, tess_epoch_ra, tess_epoch_dec)
    for i in range(len(sectors)):
        print('Sector: {:2d}, Camera: {}, CCD: {}, X pixel: {:8.3f}, Y pixel: {:8.3f}'.format(sectors[i], cameras[i], ccds[i], np.round(target_x[i],3), np.round(target_y[i],3)))
    print('')

    #Get target pixel files. 
    tpfs = tpf_grabber(tic_id, sectors)
    for i in range(len(tpfs)):
        bottom_left_x = tpfs[i].plot().get_xlim()[0]
        bottom_left_y = tpfs[i].plot().get_ylim()[0]
    
    print('')