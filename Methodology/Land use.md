## Land use indicators

Concerned impact categories:
- Land occupation, biodiversity (midpoint)
- Land transformation, biodiversity (midpoint)
- Land occupation, biodiversity (damage)
- Land transformation, biodiversity (damage)

## 1 Description of the environmental problem
Land occupation and land transformation represent major environmental pressures because they directly alter ecosystems, 
biodiversity, and the provision of natural resources and ecosystem services. Human activities such as agriculture, 
forestry, infrastructure development, and urbanization require large areas of land, often at the expense of natural 
habitats. These interventions modify ecological functions and can trigger cascading impacts on soil quality, water 
regulation, carbon storage, and species survival.

Land occupation refers to the continuous use of land for a specific purpose (e.g. cropland, pasture, urban settlements) 
without necessarily changing the underlying natural cover permanently. While land remains occupied, ecological processes 
are constrained and natural regeneration of habitats is limited, reducing the capacity of ecosystems to support native 
species.

Land transformation, on the other hand, refers to the change of land cover, such as the conversion of forests into 
agricultural fields or wetlands into urban areas. All transformations are disruptive, as it eliminates entire ecosystems,
but also lead to the development of other ecosystems. Land transformation is considered reversible in LCIA and it is
thus considered that land will thus return to its original form after a certain recovery period.

## 2 Impact pathway
The occupation or transformation of land leads to the temporary disappearance of species in a given ecosystem.

Only the impact on biodiversity is considered in concerned impact categories. In other words, the impacts of land change
leading to additional emissions of GHG is NOT covered in the present impact categories.


## 3 Technical documentation
### 3.1 Midpoint indicator(s)
The Land occupation/transformation midpoint indicators are based on the ratio of the corresponding damage indicators 
over the damage indicator of the reference value. In other words, if the value for the occupation of artificial areas
in Canada is 0.4 and that the reference value (i.e., occupation of annual crops global) is 0.7 then the midpoint will
be 0.4 / 0.7.

### 3.2 Damage indicator(s)
The local CFs per biome per land use per taxa are taken from Chaudhary et al. (2015) [http://dx.doi.org/10.1021/acs.est.5b02507].
These CFs can be found in the source database in the table "CF - regionalized - Land use - native". 
Chaudhary et al. covers 5 taxa:
- amphibians
- birds
- mammals
- plants
- reptiles <br>
and 6 land use types:
- Annual crops
- Permanent crops
- Pasture
- Forest, used, extensive
- Forest, used, intensive
- Urban areas

For each of the biome/land use type/taxa combinations, the median value is calculated.

Then, we combine the ecoregions of Olson 2001 [https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2] with
regions of ecoinvent (GIS shapefile available on their website) and land use type maps obtained from ESA World Cover 2021,
using the Google Earth engine. 
This allows us to determine an aggregated CF for a given region, per land use type. 
The aggregation between the different taxa is based on the distributions per ecoregion also provided in Chaudhary et al.
(2015). This distribution can be found in the table "SI - Land use - taxon distribution per ecoregion", while the land
use type per ecoinvent regions can be found in the table "SI - Land use - land use type per region".

The resulting national aggregated CFs can be found in table "CF - regionalized - Land use - aggregated".

For land transformation, we used the recovery times per taxa provided in Chaudhary et al. (2015). These can be found in
table "SI - Land use - recovery times". The previously obtained CFs are multiplied by these recovery times and then 
divided by two, thus yielding the CFs for land transformation.

Note that:
1. The land use type maps only provided crops in general (no differentiation between annual and permanent crops). We thus
relied on FAOSTAT data (found directly on their database which can be consulted online, but that is not citable) 
estimating the share of temporary vs permanent crops in many countries. These distributions can be found in table "SI - 
Land use - annual vs permanent crops". For countries/regions with no FOASTAT data, a global average of annual vs 
permanent crops was derived and applied.
2. An unspecified land use type is derived by averaging all land use types together, using their respective land cover areas
as weighting. The unspecified land use type is thus composed of around ~55% pasture/meadows - ~29% forests - 16% crops
(annual and permanent) - <0.01% urban.
3. A land use type "Sealed area" was created to represent urban area where all biodiversity is displaced. This is
different from the "Artificial areas" where some biodiversity persists. For instance, on a mineral extraction site,
some biodiversity persists while in an industrial zone or on the road itself, there is no biodiversity at all. Sealed 
area was thus attributed a CF of 1 PDF.
4. Whenever no data in Chaudhary et al. (2015) was available for a given land use type in a given biome (e.g.,urban 
areas in the desert), the global average value was attributed. While this might be less scientific than keeping the 
absence of value, the latter yields null CFs for certain land use types in certain countries (e.g., urban areas in 
the United Arab Emirates), which is too extreme to keep.