import galsim
import numpy as np

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.afw.math as afwMath

from astropy.stats import mad_std
from tempfile import NamedTemporaryFile

from lsst.afw.geom import makeSkyWcs
from lsst.meas.algorithms import KernelPsf
from lsst.afw.math import FixedKernel
from lsst.pipe.tasks.coaddInputRecorder import CoaddInputRecorderTask, CoaddInputRecorderConfig
from lsst.meas.algorithms import CoaddPsf, CoaddPsfConfig


def convert_to_dm(image, psf, psf_size):
    """Convert galsim image into DM format.
    Assumes the PSF is constant across the image
    Assumes the WCS is a TanWCS

    Parameters
    ----------
    image : galsim.Image
        Image to be converted
    psf : galsim.GSObject
        PSF of image
    psf_size : int
        Size of the PSF postage stamp 

    Returns
    -------
    afwImage.ExposureF
    """

    image_sizex, image_sizey = image.array.shape

    masked_image = afwImage.MaskedImageF(image_sizex, image_sizey)
    masked_image.image.array[:] = image.array

    var = mad_std(image.array)**2
    masked_image.variance.array[:] = var
    masked_image.mask.array[:] = 0
    exp = afwImage.ExposureF(masked_image)

    psf_image = galsim.ImageF(psf_size, psf_size, wcs=image.wcs)
    psf.drawImage(psf_image)
    exp_psf = KernelPsf(FixedKernel(afwImage.ImageD(psf_image.array.astype(float))))
    exp.setPsf(exp_psf)

    calib = afwImage.makePhotoCalibFromCalibZeroPoint(27)
    exp.setPhotoCalib(calib)

    # set WCS
    cd_matrix = image.wcs.cd
    crpix = geom.Point2D(image.wcs.crpix[0], image.wcs.crpix[1])
    crval = geom.SpherePoint(image.wcs.center.ra.deg, image.wcs.center.dec.deg, geom.degrees)
    wcs = makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cd_matrix)
    exp.setWcs(wcs)
    return exp


def dm_coadd(exps, global_wcs, sky_box, kernel='lanczos3'):
    """Coadd a set of exposures

    Parameters
    ----------
    exps: list of afwImage.ExposureF
        Set of exposures to coadd
    global_wcs: geom.SkyWCS
        WCS of final coadd
    sky_box: geom.BBox2I
        bounding box of the coadd
    kernel: str
        warping kernel

    Returns
    -------
    afwImage.ExposureF, list of afwImage.ExposureF
        finral image, individual warped images
    """

    # Configuration for coadd psf
    input_recorder_config = CoaddInputRecorderConfig()
    input_recorder = CoaddInputRecorderTask(config=input_recorder_config, name="dummy")
    coadd_psf_config = CoaddPsfConfig()
    coadd_psf_config.warpingKernelName = kernel

    # warping configuration
    warp_config = afwMath.Warper.ConfigClass()
    warp_config.warpingKernelName = kernel
    warper = afwMath.Warper.fromConfig(warp_config)
    wexps = []
    weight_list = []
    for i, exp in enumerate(exps):
        # Compute variance weight
        stats_ctrl = afwMath.StatisticsControl()
        stat_obj = afwMath.makeStatistics(exp.variance, exp.mask,
                                          afwMath.MEANCLIP, stats_ctrl)
        mean_var, mean_var_err = stat_obj.getResult(afwMath.MEANCLIP)
        weight = 1.0 / float(mean_var)
        weight_list.append(weight)
        wexp = warper.warpExposure(
            global_wcs, exp, maxBBox=exp.getBBox(), destBBox=sky_box)

        # Need coadd psf because psf may not be valid over the whole image
        ir_warp = input_recorder.makeCoaddTempExpRecorder(i, 1)
        good_pixel = np.sum(np.isfinite(wexp.image.array))
        ir_warp.addCalExp(exp, i, good_pixel)
        warp_psf = CoaddPsf(ir_warp.coaddInputs.ccds, global_wcs,
                            coadd_psf_config.makeControl())
        wexp.getInfo().setCoaddInputs(ir_warp.coaddInputs)
        wexp.setPsf(warp_psf)
        wexps.append(wexp)

    # combine stack images using mean
    stats_flags = afwMath.stringToStatisticsProperty("MEAN")
    stats_ctrl = afwMath.StatisticsControl()
    stats_ctrl.setWeighted(True)
    stats_ctrl.setCalcErrorFromInputVariance(True)

    masked_images = [w.getMaskedImage() for w in wexps]
    stacked_image = afwMath.statisticsStack(
        masked_images, stats_flags, stats_ctrl, weight_list, 0, 0)
    stacked_exp = afwImage.ExposureF(stacked_image)
    stacked_exp.getInfo().setCoaddInputs(input_recorder.makeCoaddInputs())
    coadd_inputs = stacked_exp.getInfo().getCoaddInputs()
    # Build coadd psf
    for wexp, weight in zip(wexps, weight_list):
        input_recorder.addVisitToCoadd(coadd_inputs, wexp, weight)
    coadd_psf = CoaddPsf(coadd_inputs.ccds, global_wcs,
                         coadd_psf_config.makeControl())
    stacked_exp.setPsf(coadd_psf)
    stacked_exp.setWcs(global_wcs)
    return stacked_exp, wexps


def load_config(config, values):
    """Load config values from a string
    This will write a temporary file and load values

    Parameters
    ----------
    config : BaseMeasurementConfig
        Configuration object to be loaded
    values : str
        Values to be set

    Returns
    -------
    None
    """
    with NamedTemporaryFile(mode='w') as tfile:
        tfile.write(values)
        tfile.flush()
        config.load(tfile.name)
