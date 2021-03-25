import copy

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack
from astropy.coordinates import match_coordinates_sky

import lsst.afw.table as afwTable
import lsst.geom as geom

from lsst.meas.base import ForcedMeasurementConfig, ForcedMeasurementTask
from lsst.dia.pipe.forcedPhotDia import DiaReferencesTask

from .utils import load_config


def create_forced(config, ref_cat, dia_exps, truth_files):
    """Run forced photometry at a set of reference positions

    This will match the forced photometry to the truth catalogs and combine them.
    The config dict has parameters relating to how to match truth to forced photometry

    Parameters
    ----------
    config : dict
        Configuration
    ref_cat : afwTable.SourceCatalog
        Reference set of objects
    dia_exps : list of afwImage.ExposureF
        Set of images to run forced photometry on
    truth_files : list of str
        Sef of truth catalogs

    Returns
    -------
    list of astropy.Table
        Matched catalog with truth information,
    list of astropy.Table
        Catalogs of truth objects with no matches
    """

    forced_config = ForcedMeasurementConfig()
    if 'forced_meas_config' in config['forced']:
        load_config(forced_config, config['forced']['forced_meas_config'])

    forced_meas = ForcedMeasurementTask(config=forced_config, refSchema=ref_cat.schema)
    ref_cat = copy.copy(ref_cat)

    if 'id_bit' in config['diff_im']:
        id_bit = int(config['diff_im']['id_bit'])
    else:
        id_bit = 16

    matches = []
    miss = []
    for i,dia_exp in enumerate(dia_exps):

        fsrc = forced_meas.generateMeasCat(dia_exp, ref_cat, dia_exp.getWcs(), afwTable.IdFactory.makeSource(i, id_bit))

        forced_meas.run(fsrc, dia_exp, ref_cat, dia_exp.getWcs())
        fcoord = SkyCoord(ra=fsrc['coord_ra']*u.rad, dec=fsrc['coord_dec']*u.rad)

        truth = Table.read(truth_files[i])
        truth = truth[truth['transient']==1]

        tx, ty = dia_exp.getWcs().skyToPixelArray(truth['ra'], truth['dec'], degrees=True)
        inside = [dia_exp.getBBox().contains(geom.Point2I(x, y)) for x,y in zip(tx, ty)]
        truth = truth[inside]

        tcoord = SkyCoord(ra=truth['ra']*u.deg, dec=truth['dec']*u.deg)
        idx, d2d, d3d = match_coordinates_sky(fcoord, tcoord)
        matched = copy.copy(fsrc).asAstropy()

        matched['d2d'] = d2d.arcsec
        matched['index'] = [i]*len(matched)
        matched['file'] = [truth_files[i]]*len(matched)
        matches.append(hstack([matched, truth[idx]]))

        tidx, td2d, td3d = match_coordinates_sky(tcoord, fcoord)
        no_match = td2d > config['diff_im']['match_dist']*u.arcsec

        missing_truth = truth[no_match]
        missing_truth['x'] = tx[inside][no_match]
        missing_truth['y'] = ty[inside][no_match]
        missing_truth['index'] = [i]*len(missing_truth)
        missing_truth['file'] = [truth_files[i]]*len(missing_truth)
        miss.append(missing_truth)

    return matches, miss

