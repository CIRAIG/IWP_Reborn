# Climate change indicators

Concerned impact categories:
- Climate change, short term (midpoint)
- Climate change, long term (midpoint)
- Climate change, human health, short term (damage)
- Climate change, human health, long term (damage)
- Climate change, ecosystem quality, short term (damage)
- Climate change, ecosystem quality, long term (damage)

For the (incl. CO2 uptake) version each of these 6 impact categories were split into 4 sub-categories:
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
sometimes also a severity factor (SF)).

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
The effect factor for Ecosystem quality comes from Thomas et al., 2004 [http://dx.doi.org/10.1038/nature02121] which provided 
an average 0.119 PDF/K effect factor. We then multiplied this value by the total surface of natural terrestrial areas 
of the world. Using the ESA GlobCover Version 2.3 Land Cover Map we obtained a surface of 2.29e13m2 which results in an 
effect factor of 2.74e12 PDF.m2.yr/K.

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
SOme LCI datasets offer the use of temporary carbon storage flows (e.g., Correction flow for delayed emission of fossil 
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
