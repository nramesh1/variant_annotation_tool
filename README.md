# Variant Annotation program

This annotation program was built using the following softwares/packages:

1. [annovar](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/)
download annovar and the following databases refGene,exac03,avsnp150,dbnsfp33a and cosmic70 from [annovar download](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)
the databases should be placed in the **annovar/humandb folder**

2. R packages
a. [dplyr](https://dplyr.tidyverse.org/)
b. [tidyr](https://tidyr.tidyverse.org/)
c. [stringr](https://stringr.tidyverse.org/)
d. [ggplot2](https://ggplot2.tidyverse.org/)
e. [patchwork](https://patchwork.data-imaginist.com/)

## Usage
The variant annotation R script is stored under the name: variant_annotation.R. 

The script works when the vcf contains two samples, defined as normal and tumor for the purposes of this task. Any variant that had multiple types and depths were split into multiple rows to have one variant type per row. These additional types maybe retained or removed from the final analysis.

The script requires the annotation R script, path to the vcf file, name of the normal sample, name of the tumor sample and path to the annovar folder. 

It can be run using the following syntax:
```
Rscript variant_annotation.R $PATH/vcf normal_sample_name tumor_sample_name $PATH/annovar
```

Example of the script usage:
```
Rscript variant_annotation.R coding_folder/Challenge_data.vcf normal vaf5 coding_folder/annovar
```
## Output files
The variant annotation script outputs a tsv file (Eg: Challenge_data_annotated.tsv) and a pdf file (Eg: Challenge_data.pdf).

The tsv file utilizes all the variants from the input vcf and annotates using annovar. It can provide a large number of annotations but they have been removed using the code for the purposes of this task. Annotations such as SIFT and Polyphen2 scores could be remove from the select statement (lines 85 - 107 in the Rscript) and analyzed for functional impact of the amino acid changes. 

The final tsv file contains the following information:

| Column Name        | Description |
| -------------------| ------------- |
| CHROM              |   |
| START              |   |
| END                |   |
| REF                |   |
| ALT                |   |
| TYPE               |   |
| DP_normal          |   |
| RO_normal          |   |
| RO_vaf5            |   |
| AO_normal          |   |
| AO_vaf5            |   |
| FREQ_RO_normal     |   |
| FREQ_RO_vaf5       |   |
| FREQ_AO_normal     |   |
| FREQ_AO_vaf5       |   |
| FREQ_AO_normal     |   |
| PERC_RO_normal     |   |
| PERC_RO_vaf5       |   |
| PERC_AO_normal     |   |
| PERC_AO_vaf5       |   |
| Func.refGene	     |   |
| Gene.refGene	     |   |
| GeneDetail.refGene |   |
| ExonicFunc.refGene |   |
| AAChange.refGene   |   |
| ExonicFunc.refGene |   |
| ExAC_ALL           |   |
| Other ExAC columns |   |
| avsnp150           |   |
| cosmic70           |   |

													