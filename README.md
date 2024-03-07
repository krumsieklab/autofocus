# AutoFocus

The repository contains scripts to replicate findings in the paper **Schweickart et al., "_AutoFocus: A hierarchical framework to explore multi-omic disease associations spanning multiple scales of biomolecular interaction"_ (2023)**

<br>

## Requirements and Installation


### Hardware Requirements
  
  The code in this repository requires only a standard computer with enough RAM to support the in-memory operations.
  
### Software Requirements
  
  This code was created with R version 4.2.2 and Rstudio Version 2023.06.2+561 and tested on macOS (Sonoma 14.2.1) with a 2.3 GHz Quad-Core Intel Core i5 CPU.
  
### Cloning the Repository from GitHub
  
  In order to install this repository as an R package, run the following command:
  
  ```
  install.packages("devtools")
  devtools::install_github(repo="krumsieklab/autofocus", subdir="autofocus")
  ```
  Alternatively, the code can be cloned using git as follows:
  
  ```
  git clone https://github.com/krumsieklab/autofocus
  ```
  
### Package Requirements
  
  The following R packages are required to run the AutoFocus scripts:


  * [cluster](https://cran.r-project.org/web/packages/cluster/index.html)
  * [dendextend](https://cran.r-project.org/web/packages/dendextend/)
  * [doParallel](https://cran.r-project.org/web/packages/doParallel/)
  * [dplyr](https://cran.r-project.org/web/packages/dplyr/)
  * [DT](https://cran.r-project.org/web/packages/DT/)
  * [foreach](https://cran.r-project.org/web/packages/foreach/)
  * [ggplot2](https://cran.r-project.org/web/packages/ggplot2/)
  * [glmnet](https://cran.r-project.org/web/packages/glmnet/)
  * [htmlwidgets](https://cran.r-project.org/web/packages/htmlwidgets/)
  * [igraph](https://cran.r-project.org/web/packages/igraph/)
  * [magrittr](https://cran.r-project.org/web/packages/magrittr/)
  * [mgm](https://cran.r-project.org/web/packages/mgm/)
  * [networkD3](https://cran.r-project.org/web/packages/networkD3/)
  * [parallel](https://www.rdocumentation.org/packages/parallel/versions/3.6.2)
  * [plotly](https://cran.r-project.org/web/packages/plotly/)
  * [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/)
  * [RhpcBLASctl](https://cran.r-project.org/web/packages/RhpcBLASctl/)
  * [reshape2](https://cran.r-project.org/web/packages/reshape2/)
  * [shiny](https://cran.r-project.org/web/packages/shiny/)
  * [shinyalert](https://cran.r-project.org/web/packages/shinyalert/)
  * [shinydashboard](https://cran.r-project.org/web/packages/shinydashboard/)
  * [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
  * [tidyverse](https://cran.r-project.org/web/packages/tidyverse/)
  

<br>

## Files in this repository

<br>

| Name | Description |
| ------ | ------ |
| R/load_xlsx.R | Loads data an annotations from an Excel file to create a SummarizedExperiment object|
| R/merge_se.R | Merges a list of SummarizedExperiment objects from different platforms perfored on the same samples|
| R/initialize.R | Initializes the AutoFocus data struct,  creating hierarchical structure, and performing association analysis of analytes with phenotype of interest |
| R/run_autofocus.R | Launches the AutoFocus Shiny app to explore results |
| R/Backend_Modules.R | Helper functions necessary for building data struct in initialize.R |
| R/Frontend_Modules.R | Helper functions necessary for visualizing results in run_autofocus.R and Shiny app |
| R/threshold_analysis.R | Performs analysis of enrichment threshold on cluster results |
| R/output_cluster_table.R | Converts clusters to table form with analyte annotations |
| inst/shiny-app/autofocus/server.R | Contains all visual elements and functionality of the Shiny app |
| inst/shiny-app/autofocus/ui.R| Shiny app UI instantiation |



<br>

## Workflow Example
  
### Data loading and concatenation
  
#### **Option 1: Through manual concatenation of data matrices** 
  AutoFocus requires three data matrices:
  
  1. Data measurements: rows represent samples, columns represent measured analytes. In the case of multiple datasets, matrices should be concatenated by binding columns and sample row order needs to be uniform.
  
  2. Sample information: Number of rows should equal the number of rows in the data measurement matrix. In the case of multiple datasets, this information should be the same for all datasets and therefore sample information from just one platform is necessary. Sample information needs to contain a sample ID column, up to two phenotypes per sample, and all desired confounders.
  
  3. Feature annotations: Number of rows should equal the number of columns in the data measurement matrix. In the case of multiple datasets, matrices should be concatenated by binding rows in the same order as was concatenated for the data measurement matrix. Feature annotations need to contain an feature ID column, and in the case of multiple datasets, a platform column is recommended.
    <br>
    <br>

R/initialize.R can now be run via "initialize_R_matrices" as follows:
    <br>
```
initialize_R_matrices(data_measurement_matrix, sample_information_matrix, feature_annotation_matrix, confounders, phenotype, ...)
```
    
 Run ?initialize_R_data_matrices in an R console for more information.

#### **Option 2 (preferred): Through loading of excel sheets and merging of SummarizedExperiments**

  The simplest way to run the AutoFocus algorithm is via datasets in individual excel sheets assuming the Excel file
contains three sheets corresponding to the assay, rowData, and colData data frames. The file must satisfy the following requirements:
  \itemize{
    \item{The data sheet must be in samples-as-columns format, with the first column of the
    data frame corresponding to the feature IDs and the first row corresponding to sample IDs.}
      \item{The feature sheet must contain at least two columns each with the feature IDs (the
       first column will be converted to rownames).}
      \item{The sample sheet must contain at least one column with sample IDs.}
 }
 <br>
 
 Create a SummarizedExperiment for each excel sheet with the following function in R/load_xlsx.R:
```
load_xls(xlsx_file_path,data_sheet_name, feature_sheet_name, sample_sheet_name)
```
 <br>
 Once this has been run for each individual platform, their respective SummarizedExperiment objects should be stored in a list and can be merged in R/merge_se.R as follows:
```
merge_se(list_of_SummarizedExperiments, sample_id_column_name, feature_id_column_name, phenotypes, confounders)
```
 Run ?merge_se in an R console for more information
 <br>
 <br>
 R/initialize.R can now be run on the resulting merged SummarizedExperiment via "initialize_R_from_SE" as follows:
```
initialize_R_matrices(merged_SummarizedExperiment, confounders, phenotype, ...)
```
  Run ?initialize_R_from_SE in an R console for more information.
 <br>
 
### Determining Enrichment Threshold
   
  The AutoFocus R struct output by either initialize_R_matrices or initialize_R_from_SE can be used to explore how enrichment threshold choice will affect the output clusters. This can be done before or after (but not while) running the Shiny app. To explore how an enrichment threeshold will affect cluster output for a specific phenotype, run
```
threshold_analysis(R_struct, phenotype)
``` 
for an output plot showing how various cluster metrics change based on threshold choice. This can be used to inform what threshold to use when running the Shiny app.
 <br>
 
### Running AutoFocus Shiny app
  The AutoFocus R struct output by either initialize_R_matrices or initialize_R_from_SE can be used to run the AutoFocus Shiny app via the "run_autofocus" function in R/run_autofocus.R. If there are feature annotations you would like to explore in your clusters, define an anno_list vector with the column names of those annotations, and then run
```
run_autofocus(R_struct, anno_list)
```  
to launch a web browser to visualize and explore results. Clusters will be determined based on a user defined enrichment threshold slider.
  <br>
  
### Table format of cluster output
  
  If you would like to explore the cluster results in data table format as opposed to in the Shiny app explorer, run the "output_cluster_table" function in R/output_cluster_table.R as follows: 
```
output_cluster_table(R_struct, threshold, phenotype, feature_name, annotation_list)
```
  If you would like all clusters for multiple phenotypes, use the phenotype = "both" option.
    

### Contact and Questions
  For issues running this code, please find contact information for the corresponding author in the original publication or on the Krumsiek Lab GitHub page.

