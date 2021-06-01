import argparse
import ee
from datetime import datetime, timezone
from task_base import HIITask


class HIIWater(HIITask):
    SCALE = 300
    GPW_CADENCE = 5
    WATER_OCCURANCE_THRESHOLD = 70
    WIDE_RIVER_MIN_WIDTH = 30  # meters
    SETTLMENT_DISTANCE_FROM_COAST = 4000  # meters
    INDIRECT_DISTANCE = 15000  # meters
    DECAY_CONSTANT = -0.0002
    INDIRECT_INFLUENCE = 4  # TODO: Set to 10.

    inputs = {
        "gsw": {
            "ee_type": HIITask.IMAGE,
            "ee_path": "JRC/GSW1_2/GlobalSurfaceWater",
            "static": True,
        },
        "caspian_sea": {
            "ee_type": HIITask.FEATURECOLLECTION,
            "ee_path": "projects/HII/v1/source/phys/caspian",
            "static": True,
        },
        "gpw": {
            "ee_type": HIITask.IMAGECOLLECTION,
            "ee_path": "CIESIN/GPWv411/GPW_Population_Density",
            "maxage": 5,  # years
        },
        "ocean": {
            "ee_type": HIITask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/ESACCI-LC-L4-WB-Ocean-Map-150m-P13Y-2000-v40",
            "static": True,
        },
        "watermask": {
            "ee_type": HIITask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "static": True,
        },
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.gpw = ee.ImageCollection(self.inputs["gpw"]["ee_path"])
        self.gsw = ee.Image(self.inputs["gsw"]["ee_path"])
        self.ocean = ee.Image(self.inputs["ocean"]["ee_path"])
        self.caspian_sea_fc = ee.FeatureCollection(
            self.inputs["caspian_sea"]["ee_path"]
        )
        self.watermask = ee.Image(self.inputs["watermask"]["ee_path"])
        self.ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
        self.kernels = {
            "coastal_settlements": ee.Kernel.euclidean(
                radius=self.SETTLMENT_DISTANCE_FROM_COAST, units="meters"
            ),
            "indirect": ee.Kernel.euclidean(
                radius=self.INDIRECT_DISTANCE, units="meters"
            ),
            "river_width_shrink": ee.Kernel.circle(
                radius=self.WIDE_RIVER_MIN_WIDTH * 2, units="meters"
            ),
            "river_width_expand": ee.Kernel.euclidean(
                radius=self.WIDE_RIVER_MIN_WIDTH * 2, units="meters"
            ),
        }

    def gpw_earliest(self):
        return self.gpw.sort("system:time_start").first()

    def gpw_latest(self):
        return self.gpw.sort("system:time_start", False).first()

    def gpw_interpolated(self):
        gpw_prev = self.gpw.filterDate(
            self.ee_taskdate.advance(-self.gpw_cadence, "year"), self.ee_taskdate
        ).first()
        gpw_next = self.gpw.filterDate(
            self.ee_taskdate, self.ee_taskdate.advance(self.gpw_cadence, "year")
        ).first()

        gpw_delta_days = gpw_next.date().difference(gpw_prev.date(), "day")
        taskdate_delta_days = self.ee_taskdate.difference(gpw_prev.date(), "day")

        gpw_diff = gpw_next.subtract(gpw_prev)

        gpw_daily_change = gpw_diff.divide(gpw_delta_days)
        gpw_change = gpw_daily_change.multiply(taskdate_delta_days)

        return gpw_prev.add(gpw_change)

    def calc(self):
        gpw_taskdate = None
        ee_taskdate_millis = self.ee_taskdate.millis()
        gpw_first_date = ee.Date(
            self.gpw.sort("system:time_start").first().get("system:time_start")
        ).millis()
        gpw_last_date = ee.Date(
            self.gpw.sort("system:time_start", False).first().get("system:time_start")
        ).millis()
        start_test = ee_taskdate_millis.lt(gpw_first_date)
        end_test = ee_taskdate_millis.gt(gpw_last_date)
        interpolate_test = start_test.eq(0).And(end_test.eq(0))
        if interpolate_test.getInfo():
            gpw_taskdate = self.gpw_interpolated()
        elif end_test.getInfo():
            gpw_taskdate = self.gpw_latest()
        elif start_test.getInfo():
            gpw_taskdate = self.gpw_earliest()
        else:
            raise Exception("no valid GPW image")

        caspian_sea = self.caspian_sea_fc.reduceToImage(["glwd_id"], ee.Reducer.max())

        ocean = ee.Image.cat(self.ocean.eq(0), caspian_sea.eq(1)).reduce(
            ee.Reducer.max()
        )

        # COASTAL
        coastal_settlements = (
            ocean.distance(self.kernels["coastal_settlements"])
            .updateMask(gpw_taskdate.gte(10))
            .selfMask()
        )

        coastal_influence = (
            coastal_settlements.distance(self.kernels["indirect"], False)
            .multiply(self.DECAY_CONSTANT)
            .exp()
            .multiply(self.INDIRECT_INFLUENCE)
            .updateMask(ocean.eq(0))
        )

        # INLAND
        surface_water = (
            self.gsw.select("occurrence")
            .lte(self.WATER_OCCURANCE_THRESHOLD)
            .unmask(1)
            .updateMask(ocean.eq(0))
        )

        wide_inland_water = (
            surface_water.reduceNeighborhood(
                reducer=ee.Reducer.max(),
                kernel=self.kernels["river_width_shrink"],
            )
            .eq(0)
            .distance(self.kernels["river_width_expand"], False)
            .gte(0)
            .selfMask()
            .reproject(crs=self.crs, scale=self.WIDE_RIVER_MIN_WIDTH)
        )

        inland_water_influence = (
            wide_inland_water.distance(self.kernels["indirect"], False)
            .multiply(self.DECAY_CONSTANT)
            .exp()
            .multiply(self.INDIRECT_INFLUENCE)
        )

        water_driver = (
            coastal_influence.unmask(0)
            .max(inland_water_influence.unmask(0))
            .selfMask()
            .updateMask(self.watermask)
            .multiply(100)
            .int()
            .rename("hii_water_driver")
        )

        self.export_image_ee(water_driver, f"driver/water")

    def check_inputs(self):
        super().check_inputs()
        # add any task-specific checks here, and set self.status = self.FAILED if any fail


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate", default=datetime.now(timezone.utc).date())
    options = parser.parse_args()
    water_task = HIIWater(**vars(options))
    water_task.run()
