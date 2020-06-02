import argparse 
import ee
from datetime import datetime, timezone
from task_base import EETask


class HIIWater(EETask):
    ee_rootdir = "projects/HII/v1/sumatra_poc"
    ee_driverdir = 'driver/water'
    # if input lives in ee, it should have an "ee_path" pointing to an ImageCollection/FeatureCollection
    inputs = {
        "jrc": {
            "ee_type": EETask.IMAGE,
            "ee_path": "JRC/GSW1_0/GlobalSurfaceWater"
        },
        "caspian": {
            "ee_type": EETask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/caspian"
        },
        "gpw": {
            "ee_type": EETask.IMAGECOLLECTION,
            "ee_path": "CIESIN/GPWv411/GPW_Population_Density",
            "maxage": 5  # years
        },
        "ocean": {
            "ee_type": EETask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/ESACCI-LC-L4-WB-Ocean-Map-150m-P13Y-2000-v40",
        },
        "watermask": {
            "ee_type": EETask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
        }
            }
    

    gpw_cadence = 5

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_aoi_from_ee("{}/sumatra_poc_aoi".format(self.ee_rootdir))


    def calc(self):

        watermask = ee.Image(self.inputs['watermask']['ee_path'])
        caspian = ee.Image(0).clip(ee.FeatureCollection(self.inputs['caspian']['ee_path'])).unmask(1).eq(0).reproject(crs='EPSG:4326',scale=300)
        ocean = ee.Image(self.inputs['ocean']['ee_path']).eq(0).add(caspian)
        jrc = ee.Image(self.inputs['jrc']['ee_path'])\
                        .select('occurrence')\
                        .lte(70)\
                        .unmask(1)\
                        .multiply(ee.Image(0).clip(ee.FeatureCollection(self.inputs['caspian']['ee_path'])).unmask(1))


        gpw = ee.ImageCollection(self.inputs['gpw']['ee_path'])
        watermask = ee.Image(self.inputs['watermask']['ee_path'])
        ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
        gpw_prior = gpw.filterDate(ee_taskdate.advance(-self.gpw_cadence, 'year'), ee_taskdate).first()
        gpw_later = gpw.filterDate(ee_taskdate, ee_taskdate.advance(self.gpw_cadence, 'year')).first()
        gpw_diff = gpw_later.subtract(gpw_prior)
        numerator = ee_taskdate.difference(gpw_prior.date(), 'day')
        gpw_diff_fraction = gpw_diff.multiply(numerator.divide(self.gpw_cadence * 365))
        gpw_taskdate = gpw_prior.add(gpw_diff_fraction)
        gpw_taskdate_300m = gpw_taskdate.resample().reproject(crs=self.crs, scale=self.scale)

        DECAY_CONSTANT = -0.0002
        INDIRECT_INFLUENCE = 4

        #COASTAL
        coast_within_4km = ocean.reduceNeighborhood(reducer=ee.Reducer.max(),kernel=ee.Kernel.circle(4000,'meters'))
        coast_within_15km = ocean.reduceNeighborhood(reducer=ee.Reducer.max(),kernel=ee.Kernel.circle(15000,'meters'))

        coastal_settlements = gpw_taskdate_300m.gte(10).multiply(coast_within_4km)

        coastal_venter = coastal_settlements.add(ee.Image(1))\
            .log()\
            .multiply(ee.Image(3.333))

        coastset_indirect = coastal_settlements.eq(0).cumulativeCost(coastal_settlements,15000).reproject(crs='EPSG:4326',scale=300).unmask(0)\
                                                .multiply(DECAY_CONSTANT).exp()\
                                                .multiply(INDIRECT_INFLUENCE)\
                                                .multiply(coast_within_15km)

        #INLAND                                            
        jrc_60_min = jrc.reduceNeighborhood(reducer=ee.Reducer.min(),kernel=ee.Kernel.circle(60,'meters')).reproject(crs='EPSG:4326',scale=30)
                                                    
        jrc_wide_rivers = jrc_60_min.reduceNeighborhood(reducer=ee.Reducer.max(),kernel=ee.Kernel.circle(60,'meters')).reproject(crs='EPSG:4326',scale=30)
                                                    
        within_4_km_of_inland_coast = jrc_wide_rivers.reduceNeighborhood(reducer=ee.Reducer.max(),kernel=ee.Kernel.circle(2000,'meters')).reproject(crs='EPSG:4326',scale=300)

        minpix_4km_threshold = within_4_km_of_inland_coast.selfMask().connectedPixelCount(50).reproject(crs='EPSG:4326',scale=900)

        inland_community_mask_4km = minpix_4km_threshold.eq(50).unmask(0).reproject(crs='EPSG:4326',scale=300)

        inland_settlements = gpw_taskdate_300m.gte(10).multiply(inland_community_mask_4km)

        inland_community_mask_11km = inland_community_mask_4km.reduceNeighborhood(reducer=ee.Reducer.max(),kernel=ee.Kernel.circle(11000,'meters'))\
                                        .reproject(crs='EPSG:4326',scale=300)

        inlandset_indirect = inland_settlements.eq(0).cumulativeCost(coastal_settlements,15000).reproject(crs='EPSG:4326',scale=300).unmask(0)\
                                        .multiply(DECAY_CONSTANT).exp()\
                                        .multiply(INDIRECT_INFLUENCE)\
                                        .multiply(inland_community_mask_11km)

        hii_water_driver = coastset_indirect.add(inlandset_indirect).updateMask(watermask).multiply(4)

        self.export_image_ee(hii_water_driver, '{}/{}'.format(self.ee_driverdir, 'hii_water_driver'))

    def check_inputs(self):
        super().check_inputs()
        # add any task-specific checks here, and set self.status = self.FAILED if any fail


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--taskdate', default=datetime.now(timezone.utc).date())
    options = parser.parse_args()
    water_task = HIIWater(**vars(options))
    water_task.run()
