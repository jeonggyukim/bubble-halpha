# save synthetic HI data to fits
from astropy.io import fits
import numpy as np
import astropy.constants as ac
import astropy.units as au
import tqdm

def calc_los_IHalpha(v_ch, j_Halpha, nH, vel,
                     ds_cgs,
                     sigma_dust=(1.0e-21*au.cm**2).cgs.value,
                     los_axis=0, verbose=True):
    """
    Function to calculate Halpha intensity along a given line of sight. Note that dust scattering is ignored and constant ionized gas temperature is assumed.

    Parameters
    ----------
    v_ch : velocity channel [km/s]
    j_Halpha : Halpha volume emissivity [erg/s/cm^3/sr]
    nH : hydrogen number density
    vel : gas velocity [km/s]
    sigma_dust : dust absorption cross section [cm^2]
    los_axis : los_axis: 0 - z, 1 - y, 2 - x

    Returns
    -------
    I : intensity
    I_thin : intensity in the optically thin limit
    """
    m_H = ac.m_p
    Tion = 8000.0*au.K
    b_kms = (np.sqrt(2.0*ac.k_B*Tion/m_H)).to('km s-1').value
    l_Halpha = (6562.8*au.angstrom).value
    c_kms = ac.c.to('km s-1').value

    I = [] # Halpha intensity
    I_thin = [] # Halpha intensity in the optically thin limit
    if verbose:
        disable = False
    else:
        disable = True

    for v_ch_ in tqdm.tqdm(v_ch, disable=disable):
        phi_v = 1/np.sqrt(np.pi)/l_Halpha*(c_kms/b_kms)*np.exp(-((v_ch_-vel)/b_kms)**2)
        dtau_los = nH*sigma_dust*ds_cgs
        # dust optical depth measured from z=zmin
        tau_cumul = dtau_los.cumsum(axis=los_axis)
        # Intensity measured by an observer at z=-\infty
        I.append(np.nansum(j_Halpha*phi_v*np.exp(-tau_cumul),
                           axis=los_axis))
        # Intensity assuming optically thin emission
        I_thin.append(np.nansum(j_Halpha*phi_v, axis=los_axis))

    I = np.array(I)
    I_thin = np.array(I_thin)

    return I,I_thin

def create_fits(domain, ytds=None, kind='pyathena'):
    hdr = fits.Header()
    if kind == 'yt':
        if ytds is None:
            raise ValueError('Parameter "kind" is set to "yt" but ytds (yt dataset) is None.')
        tMyr=ytds.current_time.to('Myr').value
        le=ytds.domain_left_edge.to('pc').value
        re=ytds.domain_right_edge.to('pc').value
        dx=(ytds.domain_width/ytds.domain_dimensions).to('pc').value
        hdr['time']=float(tMyr)
        hdr['xmin']=(le[0],'pc')
        hdr['xmax']=(re[0],'pc')
        hdr['ymin']=(le[1],'pc')
        hdr['ymax']=(re[1],'pc')
        hdr['zmin']=(le[2],'pc')
        hdr['zmax']=(re[2],'pc')
        hdr['dx']=(dx[0],'pc')
        hdr['dy']=(dx[1],'pc')
        hdr['dz']=(dx[2],'pc')
    elif kind == 'pyathena_classic':
        hdr['time']=domain['time']
        hdr['xmin']=(domain['left_edge'][0],'pc')
        hdr['xmax']=(domain['right_edge'][0],'pc')
        hdr['ymin']=(domain['left_edge'][1],'pc')
        hdr['ymax']=(domain['right_edge'][1],'pc')
        hdr['zmin']=(domain['left_edge'][2],'pc')
        hdr['zmax']=(domain['right_edge'][2],'pc')
        hdr['dx']=(domain['dx'][0],'pc')
        hdr['dy']=(domain['dx'][1],'pc')
        hdr['dz']=(domain['dx'][2],'pc')
    elif kind == 'pyathena':
        hdr['time']=domain['time']
        hdr['xmin']=(domain['le'][0],'pc')
        hdr['xmax']=(domain['re'][0],'pc')
        hdr['ymin']=(domain['le'][1],'pc')
        hdr['ymax']=(domain['re'][1],'pc')
        hdr['zmin']=(domain['le'][2],'pc')
        hdr['zmax']=(domain['re'][2],'pc')
        hdr['dx']=(domain['dx'][0],'pc')
        hdr['dy']=(domain['dx'][1],'pc')
        hdr['dz']=(domain['dx'][2],'pc')

    hdu = fits.PrimaryHDU(header=hdr)

    return hdu

def add_header_for_glue(hdu,hdr,axis='xyz'):
    for i,ax in enumerate(axis):
        hdu.header['CDELT{}'.format(i+1)]=hdr['d{}'.format(ax)]
        hdu.header['CTYPE{}'.format(i+1)]=ax
        hdu.header['CUNIT{}'.format(i+1)]=hdr.comments['d{}'.format(ax)]
        hdu.header['CRVAL{}'.format(i+1)]=hdr['{}min'.format(ax)]
        hdu.header['CRPIX{}'.format(i+1)]=hdr['{}max'.format(ax)]+hdr['{}min'.format(ax)]
    return

def save_to_fits(domain,vchannel,IHalpha,IHalpha_thin,fitsname=None):
    hdul = fits.HDUList()
    hdu = create_fits(domain)

    hdu.header['vmin']=(vchannel.min(),'km/s')
    hdu.header['vmax']=(vchannel.max(),'km/s')
    hdu.header['dv']=(vchannel[1]-vchannel[0],'km/s')

    hdul.append(hdu)
    for fdata,label in zip([IHalpha,IHalpha],
                           ['IHalpha','IHalpha_thin']):
        hdul.append(fits.ImageHDU(name=label,data=fdata))

    hdr=hdu.header
    for hdu in hdul:
        add_header_for_glue(hdu,hdr,axis='xyv')

    if fitsname is not None: hdul.writeto(fitsname,overwrite=True)
    return hdul
