![image](Report_changes/2.0/images/IW+-logo.png)

Python package creating the IMPACT World+ files directly downloadable into various software.

## What's the goal of this package?
This package is used by the IMPACT World+ team to take data from characterization models (e.g., IPCC, AWARE, Usetox,
etc.), manage and harmonize them to ultimately integrate these models in the IMPACT World+ overall methodology. All
this is done automatically and in a totally transparent/open-source manner, only relying on open-access, free data. This
package could be interesting to other users to have a precise look on what assumptions are taken by the IW+ team during
the generation of the different files for the different software platforms.

## Get started
To get started you can follow the Tutorial.ipynb jupyter notebook. You will need the source database from IW+ (SQL database)
which is downloadable here: https://doi.org/10.5281/zenodo.1488368.

## The different generated files of IW+
#### User files
After running the code (follow the Tutorial.ipynb file) you will find different versions of IW+ in the Databases folder:
- a **brightway** version (in the form of a .bw2package file), linking IW+ to the flows from biosphere3 of a selected 
brightway project. Note that this is ecoinvent-dependent, meaning that each generated IW+ files for brightway only
operates with a specific ecoivnent version. That is because biosphere3 is fixed and is ecoinvent-dependent itself. Also
note that either a bw2 or bw2.5 version of the files will be generated depending on the arguments provided by the user and
obviously of the brightway version used in the brightway project used.
- **"pure" ecoinvent** versions, linking different versions of ecoinvent to IW+. These files are available in Excel format.
- a **SimaPro** version (.csv file).
- an **openLCA** version in a .zip format (importable as a JSON-LD file in openLCA)
- an **exiobase** version, linking the environmental extensions of the exiobase GMRIO database to IW+. Either for version
3.8.2 and before and for post versions 3.9. That is because the list of environmental extensions changed after version 3.9.
- a **developer** version, this only regroups all the characterization factors of the IW+ method, in an Excel file. This 
file is mostly relevant for developers wishing to integrate IW+ in their tools.

## Mappings
Note that this package includes mappings between the flows of IW+ and the flows of the different software platforms, which
could be of interest.

## Authors
- Maxime Agez (maxime.agez@polymtl.ca)
- Elliot Muller (elliot.muller@polymtl.ca)

## Related scientific paper
Bulle, C., Margni, M., Patouillard, L. et al. IMPACT World+: a globally regionalized life cycle impact assessment 
method. Int J Life Cycle Assess 24, 1653–1674 (2019). https://doi.org/10.1007/s11367-019-01583-0