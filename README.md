## HII WATER DRIVER

## What does this task do?

This task calculates the (unitless) "influence" of navigable waterways on the terrestrial surface as one of the key drivers for a combined [Human Influence Index](https://github.com/SpeciesConservationLandscapes/task_hii_weightedsum). "Influence" is a pressure score based on proximity to a navigable waterway. Coasts, wide rivers and lakes are considered navigable if the meet key criteria related to the distance from a population center (for coastlines including that of the Caspian Sea) or based on width and connectivity (inland waters). These methods are adapted from the logic followed by [Venter et al. 2016](https://www.nature.com/articles/sdata201667)).

The source dataset for ocean coastlines and the Caspian Sea are derived from [data Adam put in... need source]. The source data for inland water is derived from [Joint Research Centre Global Surface Water (GSW)](https://global-surface-water.appspot.com/) product. The source population density cells are derived from the WoldPop Population Data dataset developed by [WorldPop](https://www.worldpop.org/). This dataset models the distribution of the global human population annually beginning in 2000 at a spatial resolution of 100 m. As a class property of HIITask the original dataset values are converted from the number of people per 100m x 100m grid cell to actual population density of people/sq km.

Inland water is defined as areas in the GSW dataset with occurrence values greater than or equal to 40. These waters are considered navigable if they have a minimum width of 30 m and are connected to at least 1024 pixels meeting the minimum width threshold. Navigable coastlines are defined as coastal areas with in 80 km from coastal settlements. Coastal settlements are any areas with population densities greater than or equal to 10 people/sq km within 4 km of a coast line.

The influence on the terrestrial surface of these navigable waterways is calculated using an exponential decay function from 0 to 15 km from the waterway. This is calculated as:

```
influence = e^(distance * decay_constant) * indirect_weight
```

## Variables and Defaults

### Environment variables

    SERVICE_ACCOUNT_KEY=<GOOGLE SERVICE ACCOUNT KEY>

### Class constants

    scale = 300
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

## Usage

    /app # python hii_popdens.py --help
    usage: task.py [-h] [-d TASKDATE]

    optional arguments:
      -h, --help            show this help message and exit
      -d TASKDATE, --taskdate TASKDATE
