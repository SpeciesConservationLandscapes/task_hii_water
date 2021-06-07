import argparse
import ee
from datetime import datetime, timezone
from task_base import HIITask


class HIIWater(HIITask):
    SETTLEMENT_DISTANCE_FROM_COAST = 4000  # meters
    COASTAL_SETTLEMENT_POPULATION_DENSITY = 10
    COASTAL_NAVIGATION_DISTANCE = 80000  # meters
    OCEAN_BUFFER_DISTANCE = 300  # meters
    GSW_OCCURRENCE_THRESHOLD = 40
    WIDE_RIVER_MIN_WIDTH = 30  # meters
    GSW_CONNECTED_PIXEL_MIN = 1024  # pixels
    INDIRECT_DISTANCE = 15000  # meters
    DECAY_CONSTANT = -0.0003
    INDIRECT_INFLUENCE = 10
    scale = 300
    gpw_cadence = 5

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
        self.gpw = ee.ImageCollection(
            self.inputs["gpw"]["ee_path"]
        )  # TODO: modify to night lights (Venter)
        self.gsw = ee.Image(self.inputs["gsw"]["ee_path"])
        self.ocean = ee.Image(self.inputs["ocean"]["ee_path"])
        self.caspian_sea_fc = ee.FeatureCollection(
            self.inputs["caspian_sea"]["ee_path"]
        )
        self.watermask = ee.Image(self.inputs["watermask"]["ee_path"])
        self.ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
        self.kernels = {
            "coastal_settlements": ee.Kernel.euclidean(
                radius=self.SETTLEMENT_DISTANCE_FROM_COAST, units="meters"
            ),
            "coastal_navigation": ee.Kernel.euclidean(
                radius=self.COASTAL_NAVIGATION_DISTANCE, units="meters"
            ),
            "ocean_buffer": ee.Kernel.euclidean(
                radius=self.OCEAN_BUFFER_DISTANCE, units="meters"
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
            .updateMask(gpw_taskdate.gte(self.COASTAL_SETTLEMENT_POPULATION_DENSITY))
            .selfMask()
        )

        coast_navigable = ocean.updateMask(
            coastal_settlements.distance(self.kernels["coastal_navigation"], False)
            .gte(0)
            .reproject(
                crs=self.crs, scale=1000
            )  # needs to be done at coarse scale to work at this distance
        )

        coast_influence = (
            coast_navigable.distance(self.kernels["indirect"], False)
            .multiply(self.DECAY_CONSTANT)
            .exp()
            .multiply(self.INDIRECT_INFLUENCE)
            .updateMask(ocean.eq(0))
        )

        # INLAND
        coast_buffer = (
            ocean.distance(self.kernels["ocean_buffer"])
            .gte(0)
            .unmask(0)
            .eq(0)
            .selfMask()
        )

        inland_water = (
            self.gsw.select("occurrence")
            .lte(self.GSW_OCCURRENCE_THRESHOLD)
            .updateMask(coast_buffer)
        )

        inland_water_mask = inland_water.eq(0).selfMask().multiply(0).unmask(1)

        inland_water_navigable = inland_water.updateMask(
            inland_water.reduceNeighborhood(
                reducer=ee.Reducer.max(),
                kernel=self.kernels["river_width_shrink"],
            )
            .eq(0)
            .distance(self.kernels["river_width_expand"], False)
            .gte(0)
        )

        inland_water_connected = (
            inland_water_navigable.connectedPixelCount(
                self.GSW_CONNECTED_PIXEL_MIN, True
            )
            .gte(self.GSW_CONNECTED_PIXEL_MIN)
            .selfMask()
            .reproject(crs=self.crs, scale=self.scale)
        )

        inland_water_influence = (
            inland_water_connected.distance(self.kernels["indirect"], False)
            .multiply(self.DECAY_CONSTANT)
            .exp()
            .multiply(self.INDIRECT_INFLUENCE)
        )

        water_driver = (
            coast_influence.unmask(0)
            .max(inland_water_influence.unmask(0))
            .selfMask()
            .updateMask(self.watermask)
            .updateMask(inland_water_mask)
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
