## HII WATER DRIVER

## What does this task do?

This task calculates the (unitless) "influence" of navigable waterways on the terrestrial surface as one of the key drivers for a combined [Human Influence Index](https://github.com/SpeciesConservationLandscapes/task_hii_weightedsum). "Influence" is a pressure score based on proximity to a navigable waterway. Coasts, wide rivers and lakes are considered navigable if the meet key criteria related to the distance from a population center (for coastlines including that of the Caspian Sea) or based on width and connectivity (inland waters). These methods are adapted from the logic followed by [Venter et al. 2016](https://www.nature.com/articles/sdata201667)).

Inland navigable waterways are based on cells in the Global Surface Water (GSW) dataset with occurrence values greater than or equal to 40. These waters are considered navigable if they have a minimum width of 30 m and are connected to at least 1024 pixels meeting the minimum width threshold. Navigable coastlines are defined as coastal areas within 80 km from coastal settlements. Coastal settlements are any areas with population densities greater than or equal to 10 people/sq km within 4 km of a coast line.

The influence on the terrestrial surface of these navigable waterways is calculated using an exponential decay function from 0 to 15 km from the waterway. This is calculated as:

```
influence = e^(distance * decay_constant) * indirect_weight
```

## Sources and preprocessing

Ocean waters are a combination of all cells with a value of 0 from the static ESA landcover water bodies data ingested into Earth Engine from  
ftp://geo10.elie.ucl.ac.be/v207/ESACCI-LC-L4-WB-Ocean-Land-Map-150m-P13Y-2000-v4.0.tif  
([Lamarche C, Santoro M, Bontemps S, D’Andrimont R, Radoux J, Giustarini L, Brockmann C, Wevers J, Defourny P, Arino O. Compilation and Validation of SAR and Optical Data Products for a Complete and Global Map of Inland/Ocean Water Tailored to the Climate Modeling Community. Remote Sensing. 2017; 9(1):36. https://doi.org/10.3390/rs9010036](https://doi.org/10.3390/rs9010036)) combined with the (static) Caspian Sea image ingested from the vector data downloaded from  
https://maps.princeton.edu/catalog/stanford-zb452vm0926  
The Caspian Sea is thus analyzed in this task in terms of coastal access, rather than inland navigation.

The source data for inland water is the [Joint Research Centre Global Surface Water (GSW)](https://global-surface-water.appspot.com/) product available within Earth Engine at `JRC/GSW1_1/GlobalSurfaceWater`. Documentation is available at [Jean-Francois Pekel, Andrew Cottam, Noel Gorelick, Alan S. Belward, High-resolution mapping of global 
surface water and its long-term changes. Nature 540, 418-422 (2016). [doi:10.1038/nature20584]](https://www.nature.com/articles/nature20584) and an accompanying [Data Users Guide](https://storage.googleapis.com/global-surface-water/downloads_ancillary/DataUsersGuidev2.pdf).

The source population density cells are derived from the WoldPop Population Data dataset developed by [WorldPop](https://www.worldpop.org/). This dataset models the distribution of the global human population annually beginning in 2000 at a spatial resolution of 100 m. As a class property of HIITask the original dataset values are converted from the number of people per 100m x 100m grid cell to actual population density of people/sq km.

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

    /app # python task.py --help
    usage: task.py [-h] [-d TASKDATE] [--overwrite]
    
    optional arguments:
      -h, --help            show this help message and exit
      -d TASKDATE, --taskdate TASKDATE
      --overwrite           overwrite existing outputs instead of incrementing
