import argparse
import ee
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
    INDIRECT_INFLUENCE = 4

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
        "ocean": {
            "ee_type": HIITask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/ESACCI-LC-L4-WB-Ocean-Map-150m-P13Y-2000-v40",
            "static": True,
        },
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.gsw = ee.Image(self.inputs["gsw"]["ee_path"])
        self.ocean = ee.Image(self.inputs["ocean"]["ee_path"])
        self.caspian_sea_fc = ee.FeatureCollection(
            self.inputs["caspian_sea"]["ee_path"]
        )
        self.kernels = {
            "coastal_settlements": ee.Kernel.euclidean(
                radius=self.SETTLEMENT_DISTANCE_FROM_COAST, units="meters"
            ),
            "coastal_navigation": ee.Kernel.euclidean(
                radius=self.COASTAL_NAVIGATION_DISTANCE, units="meters"
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

    def calc(self):
        caspian_sea = self.caspian_sea_fc.reduceToImage(["glwd_id"], ee.Reducer.max())
        ocean = ee.Image.cat(self.ocean.eq(0), caspian_sea.eq(1)).reduce(
            ee.Reducer.max()
        )

        # COASTAL
        coastal_settlements = (
            ocean.distance(self.kernels["coastal_settlements"])
            .updateMask(
                self.population_density.gte(self.COASTAL_SETTLEMENT_POPULATION_DENSITY)
            )
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
        )

        # INLAND
        inland_water = self.gsw.select("occurrence").lte(self.GSW_OCCURRENCE_THRESHOLD)

        inland_water_navigable = inland_water.updateMask(
            inland_water.reduceNeighborhood(
                reducer=ee.Reducer.max(), kernel=self.kernels["river_width_shrink"],
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
            .updateMask(self.watermask)
            .multiply(100)
            .int()
            .rename("hii_water_driver")
        )

        self.export_image_ee(water_driver, f"driver/water")

    def check_inputs(self):
        super().check_inputs()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="overwrite existing outputs instead of incrementing",
    )
    options = parser.parse_args()
    water_task = HIIWater(**vars(options))
    water_task.run()
