import argparse
import ee
from datetime import datetime, timezone
from task_base import EETask


class HIIWater(EETask):
    ee_rootdir = "projects/HII/v1/sumatra_poc"
    ee_driverdir = "driver/water"
    # TODO: figure out inputs that will be updated
    inputs = {
        "jrc": {
            "ee_type": EETask.IMAGE,
            "ee_path": "JRC/GSW1_1/GlobalSurfaceWater",
            "maxage": 30,
            "static": True,
        },
        "caspian": {
            "ee_type": EETask.FEATURECOLLECTION,
            "ee_path": "projects/HII/v1/source/phys/caspian",
            "maxage": 30,
        },
        "gpw": {
            "ee_type": EETask.IMAGECOLLECTION,
            "ee_path": f"{ee_rootdir}/misc/gpw_interpolated",
            "maxage": 1,
        },
        "ocean": {
            "ee_type": EETask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/ESACCI-LC-L4-WB-Ocean-Map-150m-P13Y-2000-v40",
            "maxage": 30,
        },
        "watermask": {
            "ee_type": EETask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "maxage": 30,
        },
    }
    scale = 300
    gpw_cadence = 5
    wide_river_width = 30  # meters
    coast_settle_min = 4000  # meters
    coast_settle_max = 15000  # meters
    DECAY_CONSTANT = -0.0002
    INDIRECT_INFLUENCE = 4

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_aoi_from_ee("{}/sumatra_poc_aoi".format(self.ee_rootdir))

    def calc(self):
        gpw, gpw_date = self.get_most_recent_image(
            ee.ImageCollection(self.inputs["gpw"]["ee_path"])
        )
        watermask = ee.Image(self.inputs["watermask"]["ee_path"])
        caspian = (
            ee.Image(0)
            .clip(ee.FeatureCollection(self.inputs["caspian"]["ee_path"]))
            .unmask(1)
        )
        ocean = (
            ee.Image(self.inputs["ocean"]["ee_path"])
            .eq(0)
            .add(caspian.eq(0).reproject(crs=self.crs, scale=self.scale))
        )
        jrc = (
            ee.Image(self.inputs["jrc"]["ee_path"])
            .select("occurrence")
            .lte(70)
            .unmask(1)
            .multiply(caspian)
        )

        # COASTAL
        coast_within_min = ocean.reduceNeighborhood(
            reducer=ee.Reducer.max(),
            kernel=ee.Kernel.circle(self.coast_settle_min, "meters"),
        )
        coast_within_max = ocean.reduceNeighborhood(
            reducer=ee.Reducer.max(),
            kernel=ee.Kernel.circle(self.coast_settle_max, "meters"),
        )
        coastal_settlements = gpw.gte(10).multiply(coast_within_min)
        # coastal_venter = (
        #     coastal_settlements.add(ee.Image(1)).log().multiply(ee.Image(3.333))
        # )

        coastset_indirect = (
            coastal_settlements.eq(0)
            .cumulativeCost(coastal_settlements, self.coast_settle_max)
            .reproject(crs=self.crs, scale=self.scale)
            .unmask(0)
            .multiply(self.DECAY_CONSTANT)
            .exp()
            .multiply(self.INDIRECT_INFLUENCE)
            .multiply(coast_within_max)
        )

        # INLAND
        jrc_min = jrc.reduceNeighborhood(
            reducer=ee.Reducer.min(),
            kernel=ee.Kernel.circle(self.wide_river_width * 2, "meters"),
        ).reproject(crs=self.crs, scale=self.wide_river_width)

        jrc_wide_rivers = jrc_min.reduceNeighborhood(
            reducer=ee.Reducer.max(),
            kernel=ee.Kernel.circle(self.wide_river_width * 2, "meters"),
        ).reproject(crs=self.crs, scale=self.wide_river_width)

        within_min_of_inland_coast = jrc_wide_rivers.reduceNeighborhood(
            reducer=ee.Reducer.max(),
            kernel=ee.Kernel.circle(self.coast_settle_max / 2, "meters"),
        ).reproject(crs=self.crs, scale=self.scale)

        # TODO: replace hardcoded values with class properties
        minpix_4km_threshold = (
            within_min_of_inland_coast.selfMask()
            .connectedPixelCount(50)
            .reproject(crs=self.crs, scale=900)
        )

        inland_community_mask_4km = (
            minpix_4km_threshold.eq(50)
            .unmask(0)
            .reproject(crs=self.crs, scale=self.scale)
        )

        inland_settlements = gpw.gte(10).multiply(inland_community_mask_4km)

        inland_community_mask = inland_community_mask_4km.reduceNeighborhood(
            reducer=ee.Reducer.max(),
            kernel=ee.Kernel.circle(
                self.coast_settle_max - self.coast_settle_min, "meters"
            ),
        ).reproject(crs=self.crs, scale=self.scale)

        inlandset_indirect = (
            inland_settlements.eq(0)
            .cumulativeCost(coastal_settlements, self.coast_settle_max)
            .reproject(crs=self.crs, scale=self.scale)
            .unmask(0)
            .multiply(self.DECAY_CONSTANT)
            .exp()
            .multiply(self.INDIRECT_INFLUENCE)
            .multiply(inland_community_mask)
        )

        # TODO: Why multiply by 4? Is this self.INDIRECT_INFLUENCE, and we should multiply only once?
        # TODO: Should we be creating watermask in this task, for use in other tasks?
        hii_water_driver = (
            coastset_indirect.add(inlandset_indirect).updateMask(watermask).multiply(4)
        )

        self.export_image_ee(
            hii_water_driver, "{}/{}".format(self.ee_driverdir, "hii_water_driver")
        )

    def check_inputs(self):
        super().check_inputs()
        # add any task-specific checks here, and set self.status = self.FAILED if any fail


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate", default=datetime.now(timezone.utc).date())
    options = parser.parse_args()
    water_task = HIIWater(**vars(options))
    water_task.run()
