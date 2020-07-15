# meso_patterns

This repository holds data, code, figures and text files for the paper "Elucidating coastal settlement patterns in Mesolithic south-eastern Norway by means of algorithmic classification and statistical modelling". The output from the GIS analyses is available in the *gis_output* folder, making it possible to skip straight to the statistical analyses in R. 

### 1. Repository content (listed in order of analysis)

    .                   
    ├── gis_data
    │   ├── sites.gpkg                # Spatial geometries for sites.
    │   ├── lakes.gpkg                # Present day lakes.
    │   ├── sediments.gpkg            # Sediment data.
    │   ├── reg_split_north.gpkg      # Polygon for the northern Gunnarsrød region.
    │   ├── reg_split_south.gpkg      # Polygon for the southern Bamble region.
    │   ├── studyarea                 # Folder holding 10 m resolution DTM tiles for study area.
    │   ├── region                    # Folder holding 25 m resolution DTM tiles for larger region.
    │   └── shoreline 
    │       ├── bamble_curve.csv      # Shoreline displacement curve for Bamble.
    │       ├── gunnarsrod_curve.csv  # Shoreline displacement curve for Gunnarsrød.
    │       ├── recode_bamble.txt     # Shoreline displacement for Bamble, formated for r.recode. 
    │       └── recode_gunnarsrod.txt # Shoreline displacement for Gunnarsrød, formated for r.recode.
    ├── python  
    │   └── analysis.py               # Python script used to run analyses in GRASS GIS.
    ├── gis_output
    │   ├── site_data.csv             # Result of GIS analysis for sites.  
    │   └── sample_data.csv           # Result of GIS analysis for random sample locations.
    │   └── maps                      # Folder with precompiled maps for visualisation.
    ├── r
    │   ├── meso_patterns.Rproj       # RStudio project. 
    │   ├── analysis.r                # Main R script.
    │   └── functions
    │       ├── functions.r           # Functions for analysis.r
    │       └── plot_functions.r      # Functions for analysis.r used to unpack and plot results.
    ├── figures                       # Folder holding all final figures used in the text. 
    └── latex                         # Folder holding raw article PDF and associated LaTeX files.

    
### 2. GIS and Python
All GIS analyses were run in GRASS GIS 7.6.1 (GRASS Development Team 2017) with the v.centerpoint extension on Windows 10. These were run using a Python 2.7 script (Python Software Foundation) with the numpy module (van der Walt et al. 2011). The coordinate system used is WGS84/UTM zone 32N (EPSG:32632). The Python script has to be run from a GRASS mapset with the correct coordinate system. The script took about 30 days to execute on a laptop computer with an Intel Core i7-8550U 1.8GHz CPU (a crash makes this estimation imprecise). Additionally, QGIS 3.6.3 (QGIS Development Team 2019) and ArcGIS 10.6 (Esri 2018) were used for some of the visualisation.

### 3. R
The statistical analyses were run using R 3.6.1 (R Core Team 2019). Additional libraries used were ggplot2 (Wickham 2016), GGally (Schloerke et al. 2018), gridExtra (Auguie 2017), cowplot (Wilke 2019), png (Urbanek 2013), car (Fox & Weisberg 2019), boot (Canty & Ripley 2019; Davison & Hinkley 1997), mice (van Buuren & Groothuis-Oudshoorn 2011), ModelMetrics (Hunt 2018), caret (Kuhn 2019) and randomForest (Liaw & Wiener 2002). The main script is *analysis.r*. This loads functions from *functions.r* and *plot_functions.r*. The complete R script ran in hours on the same laptop computer as the one used for the GIS analyses.  

### 4. Data 
Apart from changing file formats and coordinate systems, as well as reducing the extent of some data, spatial data from other sources are provided in their original form. This is with the exception of *sites.gpkg*. 

*sites.gpkg* is the site data as prepared for analysis. This is originally open data from the Directorate for Cultural Heritage which has been manually cleaned to resolve issues with a few overlapping polygons. The 'askeladden_id' column in the attributes of *sites.gpkg* can be used to extract the same original data.  
Site data: https://kartkatalog.geonorge.no/metadata/c6896f24-71f9-4203-9b6f-faf3bfe1f5ed  
Stray finds: https://kartkatalog.geonorge.no/metadata/17150d2c-b50d-4792-80f4-0cb2ec5eaa79   

The site data has additionally been given a quality score based on the information available from the database *Askeladden*, associated with the spatial geometries (see article for more on this). The database is run by the Directorate for Cultural Heritage. This is not open, but researchers, students, and cultural resource management personell can apply for access:
https://askeladden.ra.no/Askeladden (the page is currently only available in Norwegian). The correct sites can be found using the IDs provided in the column 'askeladden_id' of *sites.gpkg*.

### 5. Licenses
The following data are under [the Norwegian Licence for Open Government Data (NLOD) 2.0](https://data.norge.no/nlod/en/2.0):  
* *sites.gpkg* is slightly modified site data originally provided by the Norwegian Directorate for Cultural Heritage (see above).  
* *lakes.gpkg* is from "the Norwegian lake database" as provided by [the Norwegian Water Resources and Energy Directorate](http://nedlasting.nve.no/gis/) (page only available in Norwegian).  
* *sediments.gpkg* is from "Sediments" as provided by [the Geological Survey of Norway](https://www.ngu.no/en/topic/datasets).   
* GeoTIFF files in gis_data/studyarea are from "National elevation model DTM10" as provided by [the Norwegian Mapping Authority](https://hoydedata.no/LaserInnsyn/).

GeoTIFF files in gis_data/region are from EU-DEM v.1.1 as provided by [the European Union, Copernicus Land Monitoring Service 2018, European Environment Agency](https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1). This is open and free following the Copernicus data and information policy Regulation (EU) No 1159/2013 of 12 July 2013: https://eur-lex.europa.eu/legal-content/EN/TXT/?uri=CELEX%3A32013R1159

The following files are based on shoreline displacement data derived from Sørensen et al. 2015, available in Solheim 2017, and Sørensen et al. 2014. Both are under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/):   
* *bamble.csv*, *gunnarsrod.csv*, *recode_bamble.txt*, *recode_gunnarsrod.txt*

The remaining data in this repository are released under [CC0](http://creativecommons.org/publicdomain/zero/1.0/).  
Text under [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/).  
Code under [MIT](http://opensource.org/licenses/MIT). Year: 2020, holder: \[Omitted for blinding\]

### 6. References
Auguie, G. 2017. *gridExtra: Miscellaneous Functions for "Grid" Graphics*. R package version 2.3. https://CRAN.R-project.org/package=gridExtra 

Canty, A. & B. Ripley 2019. *boot: Bootstrap R (S-Plus) Functions*. R package version 1.3-22.

Davison, A. C. & D.V. Hinkley 1997. *Bootstrap Methods and Their Applications*. Cambridge University Press, Cambridge.

Esri 2018. ArcGIS Desktop: Release 10.6. Environmental Systems Research Institute, Redlands.

Fox, J. & S. Weisberg 2019. An R Companion to Applied Regression, Third Edition. Sage, Thousand Oaks. URL: https://socialsciences.mcmaster.ca/jfox/Books/Companion/

GRASS Development Team 2019. Geographic Resources Analysis Support System (GRASS) Software, Version 7.6. Open Source Geospatial Foundation. http://grass.osgeo.org

Hunt, T. 2018. *ModelMetrics: Rapid Calculation of Model Metrics*. R package version 1.2.2. https://CRAN.R-project.org/package=ModelMetrics

Kuhn, M. Contributions from J. Wing, S. Weston, A. Williams, C. Keefer, A. Engelhardt, T. Cooper, Z. Mayer, B. Kenkel, the R Core Team, M. Benesty, R. Lescarbeau, A. Ziem, L. Scrucca, Y. Tang, C. Candan & T. Hunt. 2019. *caret: Classification and Regression Training*. R package version 6.0-84. https://CRAN.R-project.org/package=caret

Liaw, A. & M. Wiener 2002. Classification and Regression by randomForest. *R News*, 2(3), 18--22.

Python Software Foundation. Python Language Reference, version 2.7. Available at http://www.python.org

R Core Team 2019. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna. https://www.R-project.org/

Schloerke, B, J. Crowley, D. Cook, F. Briatte, M. Marbach, E. Thoen, A. Elberg & J. Larmarange 2018. *GGally: Extension to 'ggplot2'*.  R package version 1.4.0. https://CRAN.R-project.org/package=GGally

Solheim, S. 2017. Naturvitenskap og andre ekspertanalyser. In *E18 Rugtvedt-Dørdal. Arkeologiske undersøkelser av lokaliteter fra steinalder og jernalder i Bamble kommune, Telemark fylke*, edited by S. Solheim, pp. 63-76. Portal forlag, Kristiansand. https://doi.org/10.23865/noasp.58

Sørensen, R., K. E. Henningsmoen, H. Høeg, & V. Gälman 2014. Holocene landhevningsstudier i søndre Vestfold og sørøstre Telemark – revidert kurve. In *Vestfoldbaneprosjektet. Arkeologiske undersøkelser i forbindelse med ny jernbane mellom Larvik og Porsgrunn kommune. Bind 1*, edited by S. Melvold & P. Persson, pp. 236–247. Portal forlag, Kristiansand. https://doi.org/10.23865/noasp.61

Sørensen, R., H. Høeg, & V. Gälman 2015. Revidert strandlinjeforskyvningskurve for Bamble. Technical report. Report for the Project E18 Rugtvedt-Dørdal. Museum of Cultural History, Unversity of Oslo.

Urbanek, S. 2013 *png: Read and write PNG images*. R package version 0.1-7. https://CRAN.R-project.org/package=png

van Buuren, S. & K. Groothuis-Oudshoorn 2011. mice: Multivariate Imputation by Chained Equations in R. *Journal of Statistical Software*, 45(3), 1-67. http://dx.doi.org/10.18637/jss.v045.i03

van der Walt, S., C. Colbert & G. Varoquaux 2011. The NumPy Array: A Structure for Efficient Numerical Computation, *Computing in Science & Engineering*, 13, 22-30. http://dx.doi.org/10.1109/MCSE.2011.37

Wickham, H. 2016. *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag, New York. https://ggplot2.tidyverse.org

Wilke, C. O. 2019. *cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'*. R package version 1.0.0. https://CRAN.R-project.org/package=cowplot

QGIS Development Team 2019. QGIS Geographic Information System. Open Source Geospatial Foundation Project. http://qgis.osgeo.org
