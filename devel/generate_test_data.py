import numpy as np
import os
from astropy.table import Table,vstack

# Create static and variable catalogs
np.random.seed(123457)

outdir = 'ngrid_test'
if os.path.exists(outdir) is False:
  os.makedirs(outdir)

# Number objects on a side
l = 10
n = l**2
max_ra = 20.03
min_ra = 19.97
min_dec = 0.03
max_dec = -0.03

x = np.linspace(min_ra, max_ra, l)
y = np.linspace(min_dec, max_dec, l)

ra, dec = np.meshgrid(x,y)
ra = ra.flatten()
dec = dec.flatten()

flux_min = 1e4
flux_max = 5e5

flux = [1e6]*n

data = {}
data['ra'] = ra
data['dec'] = dec
data['flux'] = flux

cat = Table(data)
cat.write(f'{outdir}/static_cat.fits',overwrite=True)


nm = 50
rac = (np.random.rand(nm)-0.5)*0.1 + (max_ra+min_ra)/2
decc = (np.random.rand(nm)-0.5)*0.1 + (max_dec+min_dec)/2

fwhm = np.random.rand(nm)*0.1 + 0.6

theta = np.deg2rad(np.random.rand(nm)*180)
pixel_scale = 0.2
mdata = {}
mdata['ra'] = rac
mdata['dec'] = decc
mdata['dudx'] = np.cos(theta) * pixel_scale
mdata['dudy'] = -np.sin(theta) * pixel_scale
mdata['dvdx'] = np.sin(theta) * pixel_scale
mdata['dvdy'] = np.cos(theta) * pixel_scale
mdata['fwhm'] = fwhm

mcat = Table(mdata)
mcat.write(f'{outdir}/static_exp.fits',overwrite=True)

galsim_config = f"""
# We again use the multiple document feature.  So start with the common information:
eval_variables :
    fimage_size : &image_size 900

#Define the PSF profile
psf :
    type : Gaussian
    sigma : {{type : Catalog , col : fwhm, num: 0}}
    flux: {{ type :
     Catalog , col : flux, num: 1}}

image :
    type: Scattered
    nobjects: {n}
    xsize: *image_size
    ysize: *image_size

    stamp_size : 64  # pixels

    wcs:
        type: Tan
        dudx : {{ type : Catalog , col : 'dudx', num: 0}}
        dudy : {{ type : Catalog , col : 'dudy', num: 0}}
        dvdx : {{ type : Catalog , col : 'dvdx', num: 0}}
        dvdy : {{ type : Catalog , col : 'dvdy', num: 0}}
        units: arcsec
        origin: center
        ra : 
          type: Degrees
          theta: {{ type : Catalog , col : 'ra', num: 0}}
        dec :
          type: Degrees
          theta: {{ type : Catalog , col : 'dec', num: 0}}


    noise :
       sky_level : 5.0e6  # ADU / arcsec^2
       index_key:  image_num

    random_seed : {{ type : Sequence , first : 553728 , index_key : 'image_num' }}
    nproc : 20

    world_pos :
        type : RADec
        ra: 
            type: Degrees
            theta: {{ type : Catalog , col : 'ra', num: 1}}
        dec:
            type: Degrees
            theta: {{ type : Catalog , col : 'dec', num: 1}}

# Define some image properties that are specific to the stamp generation phase.
stamp :
    draw_method : auto

# # Define the names and format of the output files
output :
  nimages : 20
  

input:
  index_key: image_num
  catalog: 
      - file_name: {outdir}/static_exp.fits
      - file_name: {outdir}/static_cat.fits
"""


galsim_file = open(f'{outdir}/static.yaml','w')
galsim_file.write(galsim_config)
galsim_file.close()


# Generate variable catalogs
l = 13
n = l**2
max_ra = 20.03
min_ra = 19.97
min_dec = 0.03
max_dec = -0.03


x = np.linspace(min_ra, max_ra, l)
y = np.linspace(min_dec, max_dec, l)

ra, dec = np.meshgrid(x,y)
ra = ra.flatten()
dec = dec.flatten()

static = Table.read(f'{outdir}/static_cat.fits')

flux_min = 1e4
flux_max = 5e5

flux = [1e5]*n #np.random.rand(n)*(flux_max-flux_min) + flux_min


nm = 50
for i in range(nm):
    data = {}
    data['ra'] = ra
    data['dec'] = dec
    data['transient'] = [1]*len(ra)
    data['flux'] = flux#*(1 + np.random.rand(len(flux))*0.1)


    cat = vstack([static, Table(data)])
    cat.write(f'{outdir}/var_{i}_tcat.fits',overwrite=True)



rac = (np.random.rand(nm)-0.5)*0.01 + (max_ra+min_ra)/2
decc = (np.random.rand(nm)-0.5)*0.01 + (max_dec+min_dec)/2

fwhm = np.random.rand(nm)*0.2 + 0.7

theta = np.deg2rad(np.random.rand(nm)*180)
pixel_scale = 0.2
mdata = {}
mdata['ra'] = rac
mdata['dec'] = decc
mdata['nobj'] = [len(static)+len(ra)]*len(rac)
#mdata['pixel_scale'] = [pixel_scale]*nm
mdata['dudx'] = np.cos(theta) * pixel_scale
mdata['dudy'] = -np.sin(theta) * pixel_scale
mdata['dvdx'] = np.sin(theta) * pixel_scale
mdata['dvdy'] = np.cos(theta) * pixel_scale
mdata['fwhm'] = fwhm

mcat = Table(mdata)
mcat.write(f'{outdir}/var_exp.fits',overwrite=True)

# generate galsim data
galsim_config = f"""
# We again use the multiple document feature.  So start with the common information:
eval_variables :
    fimage_size : &image_size 900

#Define the PSF profile
psf :
    type : Gaussian
    sigma : {{type : Catalog , col : fwhm, num: 0}}
    flux: {{ type : Catalog , col : flux, num: '$image_num+1'}}

image :
    type: Scattered
    nobjects: {{ type : Catalog , col : 'nobj', num: 0}}
    xsize: *image_size
    ysize: *image_size

    stamp_size : 64  # pixels

    wcs:
        type: Tan
        dudx : {{ type : Catalog , col : 'dudx', num: 0}}
        dudy : {{ type : Catalog , col : 'dudy', num: 0}}
        dvdx : {{ type : Catalog , col : 'dvdx', num: 0}}
        dvdy : {{ type : Catalog , col : 'dvdy', num: 0}}
        units: arcsec
        origin: center
        ra : 
          type: Degrees
          theta: {{ type : Catalog , col : 'ra', num: 0}}
        dec :
          type: Degrees
          theta: {{ type : Catalog , col : 'dec', num: 0}}


    noise :
       sky_level : 5.0e6  # ADU / arcsec^2
       index_key:  image_num

    random_seed : {{ type : Sequence , first : 553728 , index_key : 'image_num' }}
    nproc : 8

    world_pos :
        type : RADec
        ra: 
            type: Degrees
            #theta: {{ type : Catalog , col : 'ra', num: 1}}
            theta: {{ type : Catalog , col : 'ra', num: '$image_num + 1'}}
        dec:
            type: Degrees
            #theta: {{ type : Catalog , col : 'dec', num: 1}}
            theta: {{ type : Catalog , col : 'dec', num: '$image_num + 1'}}



# Define some image properties that are specific to the stamp generation phase.
stamp :
    draw_method : auto

# # Define the names and format of the output files
output :
  nimages : 50
  

input:
  index_key: image_num
  catalog: 
      - file_name: {outdir}/var_exp.fits
"""

galsim_file = open(f'{outdir}/var.yaml','w')
galsim_file.write(galsim_config)
galsim_file.close()