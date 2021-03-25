
import yaml
import os
import argparse

import lsst.afw.image as afwImage
from dia_sim import (create_template, create_association, create_diff_images, create_forced)

def nested_get(dic, keys):    
    "Get dictionary inside a nested dict"
    for key in keys:
        dic = dic[key]
    return dic

parser = argparse.ArgumentParser(description='Run difference image simulations')
parser.add_argument('config', default='job.yaml', help='config file')

args, unknown = parser.parse_known_args()

config = yaml.safe_load(open(args.config))

# Add commandline overides
for entry in unknown:
    key, value = entry.split('=')
    keys = key.split('.')
    d = nested_get(config, keys[:-1])
    d[keys[-1]] = value

if os.path.exists(config['global']['outdir']) is False:
    os.mkdir(config['global']['outdir'])

if 'read_file' not in config['template']:
    coadd_exp, coadd_src = create_template(config)
    coadd_exp.writeFits(os.path.join(config['global']['outdir'],
                                     config['template']['output']) + '_exp.fits')
    if config['template']['write_src']:
        coadd_src.writeFits(os.path.join(config['global']['outdir'],
                                         config['template']['output']) + '_src.fits')
else:
    coadd_exp = afwImage.ExposureF(os.path.join(config['global']['outdir'],
                                                config['template']['read_file']))
    coadd_src = None

dm_ims, dia_exps, dia_srcs, truth_files = create_diff_images(config, coadd_exp, coadd_src)

if config['diff_im']['write_dia_exp']:
    for dm_im, dia_exp,dia_src,truth in zip(dm_ims, dia_exps, dia_srcs, truth_files):
        dm_name = os.path.basename(truth).replace('.fits', '_exp.fits')
        exp_name = os.path.basename(truth).replace('.fits', '_dia_exp.fits')
        src_name = os.path.basename(truth).replace('.fits', '_dia_src.fits')

        dia_exp.writeFits(os.path.join(config['global']['outdir'], exp_name))
        dia_src.writeFits(os.path.join(config['global']['outdir'], src_name))
        dm_im.writeFits(os.path.join(config['global']['outdir'], dm_name))

ref_cat, id_catalog = create_association(config, dia_exps, dia_srcs)

if config['association']['write_ref_cat']:
    ref_cat.writeFits(os.path.join(config['global']['outdir'], 'ref_cat.fits'))

matches, misses = create_forced(config, ref_cat, dia_exps, truth_files)

if config['forced']['write_forced']:
    for match,miss,truth in zip(matches, misses, truth_files):
        match_name = os.path.basename(truth).replace('.fits', '_match.fits')
        miss_name = os.path.basename(truth).replace('.fits', '_miss.fits')

        match.write(os.path.join(config['global']['outdir'], match_name), overwrite=True)
        miss.write(os.path.join(config['global']['outdir'], miss_name), overwrite=True)
