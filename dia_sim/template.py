import yaml
import galsim

import lsst.afw.table as afwTable
import lsst.geom as geom

from lsst.meas.base import SingleFrameMeasurementConfig, SingleFrameMeasurementTask
from lsst.meas.algorithms import SourceDetectionTask, SourceDetectionConfig
from lsst.afw.geom import makeSkyWcs, makeCdMatrix

from .utils import dm_coadd, convert_to_dm, load_config

def create_template(config):
    """Create a template image from a configuration file.
    It uses the galsim config code to generate individual images

    Paramteres
    ----------
    config : 

    Returns
    -------
    afwImage.ExposureF, afwImage.SourceCatalog
        The coadded template image and the set of detected objects on the coadd
    """
    detection_config = SourceDetectionConfig()
    detection_task = SourceDetectionTask(config=detection_config)

    schema = afwTable.SourceTable.makeMinimalSchema()

    # Setup algorithms to run
    meas_config = SingleFrameMeasurementConfig()
    if 'meas_config' in config['template']:
        load_config(meas_config, config['template']['meas_config'])

    # Define slots, which we don't really use
    meas_config.slots.apFlux = None
    meas_config.slots.gaussianFlux = None
    meas_config.slots.calibFlux = None
    meas_config.slots.modelFlux = None
    meas_task = SingleFrameMeasurementTask(config=meas_config, schema=schema)
    table = afwTable.SourceTable.make(schema)

    ims = []
    dm_ims = []
    psfs = []
    srcs = []

    template_config = yaml.safe_load(open(f"{config['global']['outdir']}/{config['template']['yaml_file']}"))
    if 'nproc' in config['global']:
        template_config['image']['nproc'] = config['global']['nproc']
    for i in range(config['template']['num_images']):
        im = galsim.config.BuildImage(template_config, image_num=i)
        ims.append(im)

        psf = galsim.config.BuildGSObject(template_config, 'psf')[0]
        psfs.append(psf)
        dm_im = convert_to_dm(im, psf, config['global']['psf_stamp_size'])
        dm_ims.append(dm_im)

        det_result = detection_task.run(table, dm_im)
        src = det_result.sources
        meas_task.run(measCat=src, exposure=dm_im, exposureId=i)
        srcs.append(src)

    global_cd_matrix = makeCdMatrix(scale=config['template']['coadd_scale']*geom.arcseconds)
    crpix = geom.Point2D(config['template']['coadd_image_size']/2,
                         config['template']['coadd_image_size']/2)
    crval = geom.SpherePoint(config['template']['ra_center'],
                             config['template']['dec_center'], geom.degrees)

    global_wcs = makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=global_cd_matrix)
    sky_box = geom.Box2I(geom.Point2I(0, 0),
                         geom.Point2I(config['template']['coadd_image_size'] - 1,
                                      config['template']['coadd_image_size'] - 1))

    coadd_exp, wexps = dm_coadd(dm_ims, global_wcs, sky_box)

    det_result = detection_task.run(table, coadd_exp)
    sources = det_result.sources
    meas_task.run(measCat=sources, exposure=coadd_exp, exposureId=1000)

    return coadd_exp, sources