
import yaml
import galsim
import glob

import lsst.afw.table as afwTable
import lsst.ip.diffim as ipDiffim

from lsst.meas.base import SingleFrameMeasurementConfig, SingleFrameMeasurementTask
from lsst.meas.algorithms import SourceDetectionTask, SourceDetectionConfig

from .utils import convert_to_dm, load_config

def create_diff_images(config, template_exp, template_src):
    """Generate difference images using a template

    Parameters
    ----------
    config : dict
        Configuration
    template_exp : afwImage.ExposureF
        Template image
    template_src : afwTable.SourceCatalog
        Detections on the template image

    Returns
    -------
    list of afwImage.ExposureF
        Original images
    list of afwImage.ExposureF
        Difference images
    list of afwTable.SourceCatalog
        Catalog of detections from difference images
    list of str
        Truth file used to gnerage catlogs 
    """

    var_config = yaml.safe_load(open(f"{config['global']['outdir']}/{config['diff_im']['yaml_file']}"))

    if 'nproc' in config['global']:
        var_config['image']['nproc'] = config['global']['nproc']
    files = glob.glob(f"{config['global']['outdir']}/{config['diff_im']['files']}")
    files.sort()
    for file in files:
        var_config['input']['catalog'].append({'file_name': file})

    ims = []
    dm_ims = []
    psfs = []
    for i in range(config['diff_im']['num_images']):
        im = galsim.config.BuildImage(var_config, image_num=i)
        ims.append(im)
        psf = galsim.config.BuildGSObject(var_config, 'psf')[0]
        psfs.append(psf)
        dm_im = convert_to_dm(im, psf, config['global']['psf_stamp_size'])
        dm_ims.append(dm_im)

    detection_config = SourceDetectionConfig()
    detection_task = SourceDetectionTask(config=detection_config)

    schema = afwTable.SourceTable.makeMinimalSchema()

    # Setup algorithms to run
    meas_config = SingleFrameMeasurementConfig()
    if 'meas_config' in config['diff_im']:
        load_config(meas_config, config['diff_im']['meas_config'])
    meas_config.slots.apFlux = None
    meas_config.slots.gaussianFlux = None
    meas_config.slots.calibFlux = None
    meas_config.slots.modelFlux = None
    meas_task = SingleFrameMeasurementTask(config=meas_config, schema=schema)

    diff_config = ipDiffim.ImagePsfMatchTask.ConfigClass()

    if 'diff_im_config' in config['diff_im']:
        load_config(diff_config, config['diff_im']['diff_im_config'])
    psfm = ipDiffim.ImagePsfMatchTask(config=diff_config)

    if 'debug' in config['diff_im']:
        import lsst.log
        logger = lsst.log.Log.getDefaultLogger()
        logger.setLevel(lsst.log.TRACE)
    
    dia_exps = []
    dia_srcs = []
    truth_files = []
    for i in range(config['diff_im']['num_images']):
        try:
            
            res = psfm.subtractExposures(template_exp, dm_ims[i])
            dia_exp = res.subtractedExposure

            table = afwTable.SourceTable.make(schema)
            det_result = detection_task.run(table, dia_exp)
            dia_src = det_result.sources
            meas_task.run(measCat=dia_src, exposure=dia_exp, exposureId=i)

            dia_exps.append(dia_exp)
            dia_srcs.append(dia_src)
            truth_files.append(var_config['input']['catalog'][i + 1]['file_name'])
        except Exception as e:
            print(e)
            continue

    return dm_ims, dia_exps, dia_srcs, truth_files
