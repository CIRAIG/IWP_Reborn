## Climate change indicators

Concerned impact categories:
- Climate change, short term (midpoint)
- Climate change, long term (midpoint)
- Climate change, human health, short term (damage)
- Climate change, human health, long term (damage)
- Climate change, ecosystem quality, short term (damage)
- Climate change, ecosystem quality, long term (damage)

### 1. Midpoint indicator(s)
The _climate change, short term_ (GWP100) and _climate change, long term_ (GTP100) indicators directly come from the values 
of the IPCC AR6 report. These values can be found in the table 7.SM.7 from the Chapter 7 Supplementary Material.

A value for the substance carbon monoxide was added, based on the molar mass of carbon between CO and CO2, thus giving 
a 1.57 kgCO2eq/kgCO characterization factor.

### 2. Damage indicator(s)
The damage indicators are the product of cumulative AGTP values from each GHG and effect factors assessing the effect
of a temperature increase on either human health or ecosystem quality.

In a first step, the cumulative Absolute Global Temperature change Potential (AGTP) values are calculated for the two
time horizons (100 years for short term and 500 years for long term). To do so, we rely on the IPCC AR6 report once again.
However, the AR6 report itself does not provide enough information to calculate these values. We thus use the equations
from Thomas Gasser (gasser@iiasa.ac.at) and Yue He (heyue@iiasa.ac.at) which were the basis for the IPCC report. The
cumulative AGTP values provide a temperature increase for each kg of GHG emitted.<br>
These increases are then multiplied by effect factors either in DALY/°C or in PDF.m2.yr/°C.<br>
The effect factors for Human health come from the World Health Organization 2003 report. It includes the increased 
effects of cardiovascular diseases, diarrhea, malnutrition, floods and malaria, for a total of 1.23e7 DALY/yr/°C. <br>
The effect factor for Ecosystem quality comes from Thomas, 2004 [http://dx.doi.org/10.1038/nature02121] which resulted 
in an average 0.119 PDF/°C. We then multiplied this value by the total surface of natural terrestrial areas of the world,
which include surfaces covered by closed and open broad-leaf deciduous forest, closed and open needleleaved deciduous or
evergreen forests and mixed broadleaved and deciduous forest. Using the ESA GlobCover Version 2.3 Land Cover Map we 
obtained a surface of 2.29e13m2 which results in an effect factor of 2.74e12 PDF.m2.yr/°C
