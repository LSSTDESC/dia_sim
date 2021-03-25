
import lsst.afw.table as afwTable

from lsst.dia.pipe.simple_association import SimpleAssociationTask


def create_association(config, dia_exps, dia_srcs):
    """Create a DIAObject catalog from individual sources

    Parameters
    ----------
    config: dict
        Configuration file
    dia_exps: list of afwImage.ExposureF
        Set of images, needed for flux zeropoints
    dia_srcs: list of afwTable.SourceCatalog
        Set of individual detected sources

    Returns
    -------
    afwImage.SourceCatalog
        Associated catalog
    ParquetTable
         Table with set of ids associated with each entry
    """
    footprints = []
    for dia_src in dia_srcs:
        footprint = []
        for rec in dia_src:
            footprint.append(rec.getFootprint())
        footprints.append(footprint)

    # assume there are less than 2**bit objects for each image
    if 'id_bit' in config['diff_im']:
        id_bit = int(config['diff_im']['id_bit'])
    else:
        id_bit = 16

    # start id for diaObjects at 1000
    id_factory = afwTable.IdFactory.makeSource(1000, id_bit)
    assoc_config = SimpleAssociationTask.ConfigClass()
    assoc_config.filters = ['r']
    associator = SimpleAssociationTask(config=assoc_config)
    associator.initialize(dia_srcs[0].schema, id_factory)

    for dia_src, footprint, dia_exp in zip(dia_srcs, footprints, dia_exps):
        associator.addCatalog(dia_src, 'r', 0, 0, dia_exp.getPhotoCalib(), footprint)
    ref_cat = associator.finalize(None)
    id_catalog = associator.getObjectIds()
    return ref_cat, id_catalog
