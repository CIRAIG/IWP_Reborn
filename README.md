![image](Report_changes/2.0/images/ImpactWorld_LogoFinal.jpg)

Python package creating the IMPACT World+ files directly downloadable into various software.

## What's new with IW+?
The management of the IW+ method has been completely revamped, going from an archaic management system of Excel 
spreadsheets to a fully coded and transparent approach.

## The different files of IW+
#### Source file
IW+ relies on a Microsoft Access database to store all of its characterization factors. The database can be downloaded
here (_link to come_). This database acts as the source for all IW+ but is not linked to any LCA database/software
combinations and also does not include all relevant factors, as many are extrapolated from the few covered in the 
Access version.

That is where this package comes into play. From this Access database, the code first extrapolates all factors and then
make the links to the different LCA databases/software.

#### User files
After running the code (follow the Tutorial.ipynb file) you will find different versions of IW+ in the Databases folder:
- a brightway2 version (in the form of a .bw2package file), linking IW+ to the flows from biosphere3 of a selected 
brightway package
- "pure" ecoinvent versions, linking different versions of ecoinvent (currently 3.5 to 3.8) to IW+. These files are 
available in Excel format and as a pandas.dataframe version to easily link IW+ and ecoinvent in Python.
- a SimaPro version (.csv file), linking IW+ to **_multiple database_**. The databases covered are:
  - Agribalyse
  - Agrifootprint
  - Ecoinvent
  - ELCD
  - Industry2.0
  - US-ei2.2
  - USLCI
  - WFLDB
- an openLCA version is in the making and will be available later

#### Developer file
For entities wishing to implement IW+ in their own software/code/... a dev version is also automatically generated.
This dev version produces the full version of IW+. As said previously, the Access database does not cover all factors of
IW+. Many more are extrapolated and stored in this dev version. The dev file is available in an Excel format.

## Authors
- Maxime Agez (maxime.agez@polymtl.ca)
- Elliot Muller (elliot.muller@polymtl.ca)

## Related scientific paper
Bulle, C., Margni, M., Patouillard, L. et al. IMPACT World+: a globally regionalized life cycle impact assessment 
method. Int J Life Cycle Assess 24, 1653–1674 (2019). https://doi.org/10.1007/s11367-019-01583-0