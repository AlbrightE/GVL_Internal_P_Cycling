# GVL_Internal_P_Cycling
***Spatiotemporal variation in internal phosphorus loading, sediment characteristics, water column chemistry, and thermal mixing in a hypereutrophic reservoir in southwest Iowa, USA (2019-2020)***

**Abstract.**
The primary aim of the data product is to quantify seasonal and spatial variation in sediment phosphorus fluxes in a temperate reservoir and evaluate mechanisms responsible for instances of elevated sediment phosphorus release. We studied Green Valley Lake, a hypereutrophic reservoir in southwest Iowa, USA, from 2019 to 2020. We measured sediment phosphorus flux rates and potential explanatory variables at three sites along the longitudinal gradient of the reservoir over six sampling events during winter and summer stratification as well as mixing events in the spring, summer, and fall. Ex situ sediment core incubations were used to measure sediment P release rates under ambient temperature and dissolved oxygen conditions. Explanatory variables measured included sediment phosphorus chemistry, sediment physical characteristics, epilimnetic and hypolimnetic nutrient concentrations, and thermal stratification patterns. These data will be used to identify mechanisms driving hot spots and hot moments of sediment phosphorus release, which will contribute to our understanding of how areas of lakebed and times of the year can disproportionately influence whole-lake water chemistry.

**Creators.**
1. ***Ellen A. Albright*** UW-Madison. ORCID ID 0000-0002-6226-9158. Contact: ealbright2@wisc.edu
2. ***Dr. Grace M. Wilkinson*** UW-Madison. ORCID ID 0000-0003-4051-2249.

**Data Citation.**
Albright EA, Wilkinson GM. 2021. Spatiotemporal variation in internal phosphorus loading, sediment characteristics, water column chemistry, and thermal mixing in a hypereutrophic reservoir in southwest Iowa, USA (2019-2020) ver 1. Environmental Data Initiative. https://doi.org/10.6073/pasta/d3a70c1f0d534cca8bdebd7f7483ef38  (Accessed 2021-10-14).

**Additional Metadata.** Available in the "GVL P Cycling 2019-2020 Metadata.XML" file in the main branch of the repository. Written in Ecological Metadata Language (EML). Includes variable, unit, and code information for all data tables.

**Data Cleaning Scripts.**
1. **GVL_Incubation_CALC.R** - Calculations to determine daily sediment P flux rates and code to check data quality
      
      Input data table - GVL_Incubation_RAW.csv
      
      Output data tables - GVL_Incubation_Tidy.csv and GVL_Incubation_Sum.csv
      
2. **GVL_Sediment_CALC.R** - Calculations to determine sediment physical variables, concentrations of sediment P species, and total sediment P concentration

      Input data tables - GVL_SedimentPhys_RAW.csv, GVL_SedimentP_RAW.csv, and GVL_SedTP_RAW.csv
      
      Output data tables - GVL_SedimentPhys_Tidy.csv, GVL_SedimentP_Tidy.csv, GVL_SedimentP_Sum.csv, GVL_SedTP_Tidy.csv, and GVL_SedTP_Sum.csv
      
3. **GVL_WaterChem_CALC.R** - Code to check data quality
      
      Input data table - GVL_WaterChem_RAW.csv


**Data Analysis Scripts.**
1. **GVL_ANALYSIS_SUBMIT.R** - All analyses used to reproduce the results, tables, and figures in our manuscript "Sediment phosphorus composition controls hot spots and hot moments of internal loading in a temperate reservoir"

      Input data tables - GVL_Incubation_Sum.csv, GVL_WaterChem_RAW.csv, GVL_SedimentP_Sum.csv, GVL_SedTP_Sum.csv, GVL_bathymetry.shp, GVL_sites.csv, GVL_SedimentAreaCALC.csv
