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

The **final tsv** file contains the following information:

| Column Name        | Description                                                        |
| -------------------| ------------------------------------------------------------------ |
| CHROM              | Chromosome Number                                                  |
| START              | Start Position                                                     |
| END                | End Position                                                       |
| REF                | Reference Allele                                                   | 
| ALT                | Alternate Allele(s)                                                |
| TYPE               | Type of Variation                                                  |
| DP_normal          | Total Depth of Normal Sample                                       |
| DP_vaf5            | Total Depth of Tumor/Another Sample. vaf5 is the sample name       |
| RO_normal          | Reference Allele Depth of Normal Sample                            |
| RO_vaf5            | Reference Allele Depth of Tumor Sample                             |
| AO_normal          | Alternate Allele Depth of Normal Sample                            |
| AO_vaf5            | Alternate Allele Depth of Tumor/Another Sample                     |
| FREQ_RO_normal     | Reference Allele Frequency of Normal Sample                        |
| FREQ_RO_vaf5       | Reference Allele Frequency of Tumor/Another Sample                 |
| FREQ_AO_normal     | Alternate Allele Frequency of Normal Sample                        |
| FREQ_AO_vaf5       | Alternate Allele Frequency of Tumor/Another Sample                 |
| PERC_RO_normal     | % Reads Supporting Reference Allele of Normal Sample               |
| PERC_RO_vaf5       | % Reads Supporting Reference Allele of Tumor/Another Sample        |
| PERC_AO_normal     | % Reads Supporting Alternate Allele of Normal Sample               |
| PERC_AO_vaf5       | % Reads Supporting Alternate Allele of Tumor/Another Sample        |
| Func.refGene	     | Region where the variant is                                        |
| Gene.refGene	     | Gene name associated with the variant                              |
| GeneDetail.refGene | Additional Details such as transcript ID, protein change etc       |
| ExonicFunc.refGene | Exonic variant function                                            |
| AAChange.refGene   | Amino acid change                                                  |
| ExAC_ALL           | Allele Frequency from ExAC database                                |
| Other ExAC columns | Allele Frequency from Exac database for several population groups  |
| avsnp150           | dbSNP150 annotation                                                |
| cosmic70           | COSMIC gene annotation                                             |

The **final pdf** file contains examples of visualizing the distributions and trends of some of the variables listed in the above table.

													