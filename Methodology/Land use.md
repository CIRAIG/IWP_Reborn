## Land use indicators

Concerned impact categories:
- Land occupation, biodiversity (midpoint)
- Land transformation, biodiversity (midpoint)
- Land occupation, biodiversity (damage)
- Land transformation, biodiversity (damage)


### 1. Midpoint indicator(s)
The Land occupation/transformation midpoint indicators are based on the ratio of the corresponding damage indicators 
over the damage indicator of the reference value. In other words, if the value for the occupation of artificial areas
in Canada is 0.4 and that the reference value (i.e., occupation of annual crops global) is 0.7 then the midpoint will
be 0.4 / 0.7.

### 2. Damage indicator(s)
The Land occupation, biodiversity damage indicator values are based on "de Baan, Alkemade et al." (2013) Table 1 of the 
supplementary information [https://doi.org/10.1007/s11367-012-0412-0]. This describes the impacts per type of 
occupation of land per biome (e.g., tropical, temperate, desert, etc.). These values are then linked to the grid model
of Olson et al. (2001) [https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2], where the information on 
countries, biomes and types of land occupied is intersected. For each country we are thus able to determine the 
average occupation CF for each type of land occupied, based on the biome composition of each country.

The Land transformation, biodiversity damage indicator is based on the Land occupation, biodiversity damage indicator,
multiplied by a recovery time. These recovery times are biome-specific and come from "de Baan, Mutel et al." (2013) 
Table S3 of the supplementary information [https://doi.org/10.1021/es400592q]. The biome occupation damage indicators 
are then multiplied by the corresponding biome recovery time and divided by 2.