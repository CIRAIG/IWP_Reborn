# Climate change indicators

Concerned impact categories:
- Climate change, short term (midpoint)
- Climate change, long term (midpoint)
- Climate change, human health, short term (damage)
- Climate change, human health, long term (damage)
- Climate change, ecosystem quality, marine ecosystem, short term (damage)
- Climate change, ecosystem quality, marine ecosystem,long term (damage)
- Climate change, ecosystem quality, terrestrial ecosystem, short term (damage)
- Climate change, ecosystem quality, terrestrial ecosystem,long term (damage)

For the (incl. CO2 uptake) version each of these 8 impact categories were split into 4 sub-categories:
- biogenic
- CO2 uptake
- fossil
- land transformation

## 1 Description of the environmental problem
Climate change receives a significant attention throughout the world because of the numerous and far-reaching 
consequences of a warming climate system on human health, ecosystems, economy and societies. Anthropogenic global 
warming is mainly due to the emission of greenhouse gases (GHG) from human activities. Other climate forcing agents 
such as near-term climate forcers (e.g. black carbon, NOx) or changes in land covers (terrestrial albedo) caused by 
human activities are not currently covered in IW+ v2.2, but will be in the future.

GHGs absorb infrared radiations emitted by the Earth’s surface, keeping thermal energy in the low atmosphere. This 
greenhouse effect is responsible for the Earth’s temperature regulation. Atmospheric GHG concentration is increasing 
since the last century, mainly due to the widespread use of fossil fuels, which results in the rising of the average 
temperature of the Earth’s atmosphere and oceans. The most common GHGs are carbon dioxide (CO2), methane (CH4) and 
nitrous oxide (N2O). Several ozone-depleting substances such as CFCs are also greenhouse gases, as well as different 
fluorinated gases.

The warming of the climate system leads to different disturbances in oceanic circulation schemes, extent of glaciers, 
importance of precipitations, sea level, etc. that, in turn, lead to a large number of consequences such as floods, 
desertification, diseases, extreme weather events, etc.

## 2 Impact pathway
A GHG emission first leads to an increase in atmospheric concentration, which then leads to an increase in radiative 
forcing. The GWP values used as mid-point characterization factors are at the radiative forcing level of the 
cause-effect chain.

This increase in radiative forcing causes different climate effects such as an increase in temperature, 
precipitation changes, etc. We modeled only the temperature increase pathway, the other types of climate effects were 
ignored due to lack of data and knowledge.

An increase in the Earth’s average temperature leads to potential impacts on humans and ecosystems through several pathways.

For the impact of climate change on human health, we considered 5 pathways, that is, the life years in good health lost
from the following 6 effects, whose effect grows more intense because of global warming:
- Dengue
- Diarrhoeal disease
- Malaria
- Undernutrition
- Coastal floods
- Heat/cold stress

For the impact of climate change on ecosystem quality, we considered a single pathway, that is, the temporary 
disappearance of species from a given location due to temperature increases.

## 3 Technical documentation
### 3.1 Midpoint indicator
The _climate change, short term_ (GWP100) and _climate change, long term_ (GTP100) indicators directly come from the values 
of the IPCC AR6 report. These values can be found in the table 7.SM.7 from the Chapter 7 Supplementary Material.
All these values can be found in table "CF - not regionalized - ClimateChange" of the source database.

A value for the substance carbon monoxide was added, based on the molecular mass of carbon between CO and CO2, thus 
giving a 1.57 kgCO2eq/kgCO characterization factor.

### 3.2 Damage indicator(s)
As a reminder, a characterization factor (CF) is the multiplication of a Fate factor (FF) and an Effect factor (EF) (and 
sometimes also a severity factor (SF) and an exposition factor (XF)).

#### 3.2.1 Fate factor
The Fate factor is identical for both damage indicators (i.e., for both human health and ecosystem quality). It
corresponds to the cumulative AGTP (Absolute Global Temperature change Potential) values from each GHG. An AGTP value
at a given time corresponds to the absolute change in global temperature resulting from a pulse emission of a specific 
GHG. It is expressed in K/kg GHG (more like in pK/kg given the often very low values). Cumulative AGTPs are simply the
addition of all AGTPs in a given time period. They are thus expressed in K.yr/kg.

Hence, for our short term damage impact indicators, this time period is considered from 0 to 100 years. For the long
term damage indicators, the time period is considered to be from 100 to 500 years.

So how are these AGTP values calculated?

The IPCC AR6 report Chapter 7 only provides the direct AGTP values at specific years (AGTP-50 and AGTP-100) and does not
directly provide equations required to calculated AGTP values at any given year. Thankfully, there is more information
on the official IPCC Github (https://github.com/IPCC-WG1/Chapter-7). More specifically, in the src/ar6/metrics folder,
we have access to the equations to calculate the AGTPs at any year for CO2, CH4, N2O and any other GHG (halogen_generic).
These equations were implemented in the code of IMPACT World+, inside the ```load_climate_change_cfs()``` function. To obtain the 
same meinshausen function as used in the generation of the IPCC AGTP values, we relied on the v1.6.2 of the fair Python 
package ```pip install fair==1.6.2```.

From these equations we can recalculate the AGTP-50 and AGTP-100 values provided in the AR6 Chapter-7 (<0.2% difference).

The calculated AGTP values at each year (AGTP1, AGTP2, etc. until AGTP500) and for each GHG are provided in the source
database of IW+ in the table SI - Climate change - fate factors (K/kg). Note that some GHGs had radiative efficiency values
and/or lifetime values so low that it resulted in null AGTP values. For instance, for butane the radiative efficiency is smaller
than 0.001 W/m2/ppb and has thus been rounded to 0 in the IPCC AR6 report, thus rendering impossible the calculation of
AGTP values for this compound.

#### 3.2.2 Effect factors
Fate factors in K.yr/kg GHG are to be multiplied with effect factors. These effect factors differ between both damage 
indicators. For human health it is expressed in DALY/(K.yr) and for ecosystem quality it is expressed in PDF.m2.yr/(K.yr).

##### 3.2.2.1 Human health
The effect factor for human health comes from two main sources: 1. the World Health Organization's 2014 report 
[https://iris.who.int/bitstream/handle/10665/134014/9789241507691_eng.pdf], which includes the increased risk of death 
from dengue, diarrhea, malnutrition, coastal flood and malaria, and 2. Rupcic, L. (2022) 
[https://backend.orbit.dtu.dk/ws/portalfiles/portal/329521472/PhD_Thesis_Lea_Rupcic.pdf], which includes the increased 
risk of death from heat stress.

These increased risks of death were then multiplied by the burden per death for each cause in the investigated regions and 
populations (calculated from the Global Burden of Disease Study 2021, IHME, 2024 [https://vizhub.healthdata.org/gbd-results/]).
Finally, we calculated the total DALY/(K.yr) by summing the values across all regions and causes of death, resulting in 
a total of 3.69e7 DALY/(K.yr). The weight of each risk and the exact total value can be seen in the 
"SI - Climate change - effect factors" table from the source database.

##### 3.2.2.2 Ecosystem quality
The effect factor for Ecosystem quality is based on Potentially Affected Fraction of species (PAF). A species is 
considered affected whenever their niche habitat temperature is exceeded for 5 consecutive years. This definition is
based on the article of Trisos et al. 2020 (https://doi.org/10.1038/s41586-020-2189-9) which is also used by 
Iordan-Vasquez et al. 2023 (https://doi.org/10.1016/j.resconrec.2023.107159). The latter is the basis of the equivalent
climate change ecosystem quality indicator in the GLAM methodology.

##### 3.2.2.2.1 Data sources
The niche temperature of each species is defined "using the	maximum	mean annual	temperature experienced	across the 
entire geographic range of a species, between 1850 and 2005" (Trisos et al. 2020). In other words, the authors had 155 mean annual 
temperatures (at grid level across the globe) and derived a maximum out of these 155 values. Crossing these local
temperatures with the locations where each species resides, the niche temperature of each species was defined. Same 
species living in different places thus share the same niche temperature. 
Trisos 2020 provides these temperatures for 30,652 species across 9 taxa. 4 taxa living on land: amphibians, birds, 
mammals, reptiles and 5 taxa living in the ocean: benthic marine invertebrates, corals & seagrasses, marine fishes, 
marine mammals and marine reptiles.

The distribution of species (where each species resides) is taken from the data of Iordan-Vasquez et al. (2023). The
data is curated as 100km x 100km grid cells. This covers 26,648 species: birds (n = 7177), terrestrial mammals (n = 5160), 
terrestrial reptiles (n = 4599), amphibians (n = 5998), marine mammals (n = 117), marine reptiles (n = 61), 
marine fish (n = 1822), benthic marine invertebrates (n = 916), and corals and seagrasses (n = 798). The data originally
comes from the International Union for Conservation of Nature (IUCN) and BirdLife International.

Since PAFs are defined as the exceedance of the niche temperature limit of each species over 5 consecutive years, we
need climate models to determine what the temperature will be over the next years. To do so, we re-use the data from
Iordan-Vasquez et al. (2023). They averaged 5 climate models (CESM1(CAM5), HadGEM2-ES, IPSL-CM5A-MR, MIROC5, MPI-ESM-MR)
to obtain temperatures from 2010 to 2100, at a 1.875° resolution, for three climate scenarios RCP2.6, 4.5 and 8.5. This
provides temperatures for air and ocean surface. The latter is used as a proxy for the total ocean temperature.

Finally, since we need to distinguish between land and marine species, we need to distinguish between land and ocean areas.
The maps to do so come from using MODIS MCD12Q1 land water mask through the Google Earth Engine.

##### 3.2.2.2.2 Methodology for calculating effect factors
To determine the effect factors, we first calculate the PAF for each grid cell, that is, the number of affected species
per grid cell, over the total number of species in that grid cell. In each grid cell, we thus compared the maximum estimated
temperature over the time period 2010-2100. This maximum often occurs at the 2100 year, but not always. This maximum of 
temperature is then compared to the niche temperatures of all species in the grid cell, and species for which this 
temperature is exceeded are considered affected. 
Note that since the temperature estimates at provided by increment of 
average of 5 years (e.g., the average temperature from 2050 and 2055), which means that as soon as the estimated temperature
exceeds the niche temperature, it directly corresponds to our definition of affected species (i.e., temperature exceeded 
for5 consecutive years).
Also note that for terrestrial species, the air temperature estimates are used while the estimated surface ocean temperatures
are used for marine species.

Once we have our PAFs per grid cell, we multiply them by their corresponding surface area (either land or water) in
each grid cell.

We finally divide these PAF.m2 by the temperature increase in the corresponding cell, that is, the difference between the
maximum estimated temperature over the period 2010-2100 and the temperature at 2010 (the reference year). This provides
PAF.m2/K effect factors.

At this point, we have what we call average effect factors, for each of the three RCP for which temperature
estimates are provided (RCP2.6, 4.5 and 8.5). In IMPACT World+, characterization factors are typically marginal 
characterization factors. So, we need a few more steps to get there.

Marginal effect factors are typically defined to represent the effect of a marginal increase of a pollutant in
the system. For this indicator, the marginal increase of pollutant would be an increase in emissions of GHGs, which can
be translated as an increase of temperature in our system. Therefore, to derive marginal effect factors, we
use the previously obtained average effect factors as a reference state, and proceed to add marginal increments
of temperature to the system to see what is the effect that will result from these increments.

Concretely, we follow the same methodology, except we compare the niche temperature to the estimated air/surface
ocean temperature from the climate models + an increment of temperature increase (say +0.01K). This increment must be
distributed across the different grid cells across the global, because simply adding the increment equally in each grid
cell does not account for the fact that some areas in the world are heating up faster than others. So we determine the
temperature increase between 2100 and 2010 in each grid cell. This provides us a distribution, reflecting some areas
heat up faster than others. This distribution is then applied to the temperature increment to distribute it in a logical
manner across the globe.
Similarly to the average effect factor, we then obtain PAFs, that we convert to PAF.m2 and PAF.m2/K. Finally,
we calculate the difference between the effect factor obtained through this increment and its reference state (the 
average effect factor). The obtained value is then divided by the temperature increment, which yields the final effect factor.

Since we both cover terrestrial and marine species, we decided to split the climate change, ecosystem quality damage
indicator into two indicators:
- climate change, ecosystem quality, terrestrial ecosystem
- climate change, ecosystem quality, marine ecosystem

In the end, the respective effect factors are 4.35e12 PDF.m2/K for terrestrial species and 31.3e12 PDF.m2/K for marine
species. Marine species are thus dramatically more affected than terrestrial species by climate change, which is an 
observation that can be found in other scientific articles (https://doi.org/10.1038/s41586-019-1132-4). For reference, 
the previous effect factors in IW+ v2.1 was 2.73e12 PDF.m2/K. This factor only represented terrestrial species. The 
effect of climate change on ecosystem quality overall has thus gone from 2.73e12 to 35.65e12, that is a 13-fold increase!

##### 3.2.2.2.3 PAF to PDF
Here we consider that any affected species will temporarily disappear. We thus have a 1:1 ratio conversion between PAF
and PDF. While this might be surprising, the original authors Trisos et al. (2020) state having tested how would their
own results be affected by taking a 20 consecutive years period (instead of 5) and noted that it did not significantly
impact results.

##### 3.2.2.2.4 The choice of the RCP
The calculations for the ecosystem quality effect factor were performed for the three RCPs for which the climate models
were derived by Iordan-Vasquez et al. (2023). In fine, we chose the RCP4.5 to determine the effect factors. The RCP2.6
corresponds to a very optimistic vision of the world, while the RCP8.5 could be considered a pessimistic vision. Note
that the arithmetic average of these three RCPs yielded temperatures estimates extremely close to the RCP4.5 estimates.
Also note that the effect factors obtained from the RCP4.5 are higher than the ones obtained from the RCP8.5. Choosing
the RCP4.5 thus corresponds to a conservative approach, where we take the worst case.

##### 3.2.2.2.5 The choice of the temperature increment
We selected the +0.01K temperature increment for the calculation of our marginal effect factor. We tested temperature
increments of +0.1K and +1K, which yielded similar effect factors (difference of less than 25%).

### 3.3 The case of biogenic carbon
With the v2.1 update, users have the possibility to choose between two versions of IMPACT World+ corresponding to two 
different philosophies when it comes to quantifying the effects of biogenic carbon. 

The carbon neutrality approach, which has always been the default option in IMPACT World+, characterizes the effect of 
biogenic CO2 as neutral. Hence, for GWP100, CO2 biogenic flows have a null characterization factor while biogenic CH4 has a 
lower CF than the fossil one. Considering the direct differentiation of CO into CO2, biogenic CO also has a null 
characterization factor.

In the (incl. CO2 uptake) files, biogenic carbon is accounted just as fossil carbon. However, the uptake of carbon is 
characterized with a negative CF (-1) while the release is characterized with a positive CF (+1). Other biogenic GHGs 
(CH4 and CO) are also accounted as fossil GHGs, with the same characterization factor.

Within the (incl. CO2 uptake) version, each impact category related to climate change was separated into four 
sub-categories: fossil, biogenic, land transformation, CO2 uptake. This separation is recommended by various standards 
such as the PEF and GHG protocol.

### 3.4 The case of temporary storage of carbon
Some LCI datasets offer the use of temporary carbon storage flows (e.g., Correction flow for delayed emission of fossil 
carbon dioxide). These flows are expressed in kgy. They are only characterized at the damage level (so for human health 
and ecosystem quality). They are attributed a negative CF for their short term impact, equal to the CF of their 
corresponding flow divided by 100 (because the flows are per year and our climate change short term impact is on 100 
years). So for CO2, the human health short term damage CF os 1.59e-6, the CF for temporary storage of CO2 flows will be
-1.59e-8. Then, for each year it is stored, the short term impact will be reduced. However, their long term impact CF
is equal to the exact same amount, but with a reversed sign. So it would be +1.59e-8 in our previous example. Hence, 
the temporary storage of GHGs, if looking at both short and long term impact categories, will result in a null impact,
as it is considered that this temporary delay is not helping to reduce climate change over 500 years. However, if only 
looking at the short term impact, then temporary storage has a beneficial effect. Note that the IMPACT World+ team 
recommends to always look at both short and long term impacts, for all categories.
