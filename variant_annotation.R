args = commandArgs(trailingOnly=TRUE)

#required packages
packages = c('readr', 'dplyr', 'tidyr', 'stringr', 'ggplot2', 'patchwork')
BiocManager::install(setdiff(packages, rownames(installed.packages())))
lapply(packages, require, character.only = TRUE)

## arguments to be input 
input_vcf = args[1]
normal_sample_name = args[2]
tumor_sample_name = args[3]
annovar_path = args[4]
                          
## read vcf, rename columns, convert read info columns into one column, extracted necessary info
## and convert the read info into two columns
vcf = readr::read_delim(file = input_vcf, delim = '\t', comment = '##', 
                        progress = show_progress(), 
                        col_types = cols(`#CHROM` = col_character(),
                                         POS = col_double(),
                                         ID = col_character(),
                                         REF = col_character(),
                                         ALT = col_character(),
                                         QUAL = col_double(),
                                         FILTER = col_character(),
                                         INFO = col_character(),
                                         FORMAT = col_character(),
                                         .default = col_character())) %>%
  dplyr::rename(CHROM = `#CHROM`) %>%
  tidyr::pivot_longer(cols = c(all_of(normal_sample_name), all_of(tumor_sample_name)),
                      names_to = 'SAMPLES',
                      values_to = 'SAMPLE_READ_INFO') %>%
  tidyr::separate(col = SAMPLE_READ_INFO, into = c('GT','GQ','DP','DPR','RO','QR','AO','QA'), sep = ':') %>%
  tidyr::separate(col = INFO, into = c('INFO', 'TYPE'), sep = 'TYPE=', remove = FALSE) %>%
  dplyr::select(-ID, -QUAL, -FILTER, -INFO, -FORMAT, -GT, -GQ, -DPR, -QR, -QA) %>%
  tidyr::separate_rows(data = ., ALT, TYPE, AO, sep = ',', convert = TRUE) %>%
  dplyr::mutate(RO = as.numeric(RO),
                AO = as.numeric(AO),
                FREQ_RO = as.numeric(RO/(RO+AO)),
                FREQ_AO = as.numeric(AO/(RO+AO)),
                PERC_RO = as.numeric(FREQ_RO*100), 
                PERC_AO = as.numeric(FREQ_AO*100)) %>%
  tidyr::pivot_wider(names_from = 'SAMPLES',
                     values_from = c('DP', 'RO', 'AO', 
                                     'FREQ_RO', 'FREQ_AO', 
                                     'PERC_RO', 'PERC_AO'))

#extracting file name and file path
file_name = stringr::str_split(string = basename(input_vcf), pattern = '.vcf')[[1]][1]
file_path = dirname(input_vcf)

message(paste0('reading vcf file:', input_vcf, ' and setting file names to ', file_name))

#creating annovar folder/checking if it exists
annovar_output = paste0(file_path, '/annovar_output/')
if (!dir.exists(annovar_output)){
  message('note: creating folder to save annovar files')
  dir.create(annovar_output)} else {message('note: folder to save annovar files already exists!')}

message('.....running annovar.....')
convert2annovar = paste0(annovar_path, '/convert2annovar.pl')
cmd_convert2annovar = paste0(convert2annovar, ' --format vcf4old ', input_vcf, ' -outfile ', annovar_output, file_name, '.annovar -withzyg')
system(command = cmd_convert2annovar, intern = TRUE, timeout = 0)

annotate_variation = paste0(annovar_path, '/annotate_variation.pl')
cmd_annotate_variation = paste0(annotate_variation, ' -geneanno -buildver hg19 ', annovar_output, file_name, '.annovar ', paste0(annovar_path, '/humandb/'))
system(command = cmd_annotate_variation, intern = TRUE, timeout = 0)

table_annovar = paste0(annovar_path, '/table_annovar.pl')
cmd_table_annovar = paste0(table_annovar, ' ', annovar_output, file_name, '.annovar ', paste0(annovar_path, '/humandb/ -buildver hg19 -out ', annovar_output, file_name, ' -remove -protocol refGene,exac03,avsnp150,dbnsfp33a,cosmic70 -operation g,f,f,f,f -nastring . -polish'))
system(command = cmd_table_annovar, intern = TRUE, timeout = 0)

## read annovar output file
annovar_multianno = read.delim(file = paste0(annovar_output, file_name, '.hg19_multianno.txt'), header = TRUE, sep = '\t')

## merge the vcf file in table form with annovar output
final_annotation_table = dplyr::left_join(x = vcf, y = annovar_multianno,
                                          by = c('CHROM' = 'Chr', 'POS' = 'Start',
                                                 'REF' = 'Ref', 'ALT' = 'Alt')) %>%
  dplyr::rename(START = POS) %>%
  dplyr::mutate(END = End,
                Func.refGene = factor(Func.refGene, levels = levels(addNA(Func.refGene)), 
                                      labels = c(levels(Func.refGene), 'unknown'), 
                                      exclude = NULL)) %>%
  dplyr::relocate(END, .after = START)

## output the vcf + annovar merged file
readr::write_delim(x = final_annotation_table, path = paste0(file_path, '/', file_name, '_annotated.tsv'), quote_escape = 'double', delim = '\t')

## making plots using ggplot
message('creating a few plots using ggplot')

## setting theme
my_theme = theme(axis.ticks.length = unit(5,'pt'),
                 axis.text.x = element_text(angle = 0,
                                            vjust = .5,
                                            size = 12),
                 axis.text.y = element_text(size = 12),
                 axis.title = element_text(size = 20),
                 plot.title = element_text(size = 20),
                 legend.text = element_text(size = 12),
                 legend.title = element_text(size = 12),
                 strip.text.y = element_blank(), #remove side boxes with labels
                 panel.border = element_rect(colour = 'black', fill=NA, size=1),
                 panel.spacing = unit(0, 'lines'),
                 legend.position = 'right',
                 legend.direction = 'vertical',
                 plot.margin = margin(0.1,0.1,0.1,0.1,'in'))

## distribution of number of variants by chromosome
plot_variants_summary_box = final_annotation_table %>%
  dplyr::select(CHROM, TYPE) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(VARIANTS = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(VARIANTS)) %>%
  ggplot(mapping = aes(x = '', y = VARIANTS)) +
  geom_boxplot(position = 'dodge', fill = 'cyan', outlier.shape = NA) + 
  geom_jitter(color = 'red', width = 0.1, size = 5, alpha = 1) +
  labs(title = 'distribution of number of variants by chromosome',
       x = 'chromosomes', y = 'number of variants') +
  theme_classic() + 
  my_theme

# distribution of type of variant
plot_type_col = final_annotation_table %>%
  dplyr::select(TYPE) %>%
  dplyr::group_by(TYPE) %>%
  dplyr::summarise(VARIANTS = n()) %>%
  dplyr::arrange(VARIANTS) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(VARIANTS)) %>%
  dplyr::mutate(TYPE = factor(TYPE, levels = TYPE)) %>%
  ggplot(mapping = aes(x = TYPE, y = VARIANTS)) +
  geom_col(position = 'dodge', fill = '#6495ed') + 
  labs(title = 'distribution of type of variant',
       x = 'type of variant', y = 'number of variants') +
  theme_classic() + 
  my_theme + 
  theme(axis.text.x = element_text(angle = 45))

## distribution of variants by chromosome
plot_variants_summary_col = final_annotation_table %>%
  dplyr::select(CHROM, TYPE) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(VARIANTS = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(VARIANTS)) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = CHROM)) %>%
  ggplot(mapping = aes(x = CHROM, y = VARIANTS)) +
  geom_col(position = 'dodge', fill = '#DC143C') + 
  facet_grid(~CHROM, scales = 'free_x', space = 'free_x') +
  labs(title = 'distribution of variants by chromosome',
       x = 'chromosomes', y = 'number of variants') +
  theme_classic() + 
  my_theme +
  theme(strip.background = element_rect(fill = 'cornsilk2'),
        strip.text.x = element_text(size = 12))

# distribution of effect of variant
plot_funcrefgene_col = final_annotation_table %>%
  dplyr::select(Func.refGene) %>%
  dplyr::group_by(Func.refGene) %>%
  dplyr::summarise(VARIANTS = n()) %>%
  dplyr::arrange(VARIANTS) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Func.refGene = factor(Func.refGene, levels = Func.refGene)) %>%
  ggplot(mapping = aes(x = Func.refGene, y = VARIANTS)) +
  geom_col(position = 'dodge', fill = '#6495ed') + 
  labs(title = 'distribution of effect of variant',
       x = 'effect of variant', y = 'number of variants') +
  theme_classic() + 
  my_theme + 
  theme(axis.text.x = element_text(angle = 45))

#Q1 = ((plot_variants_summary_box | plot_type_col) / (plot_variants_summary_col | plot_funcrefgene_col)) + 
#  patchwork::plot_annotation(title = 'Q1')

#depth of sequence coverage at the site of variation
plot_depth_dist = final_annotation_table %>%
  dplyr::select(paste0('DP_', normal_sample_name), paste0('DP_', tumor_sample_name)) %>%
  tidyr::pivot_longer(cols = c(paste0('DP_', normal_sample_name), paste0('DP_', tumor_sample_name)),
                      names_to = 'SAMPLE', values_to = 'DP') %>%
  dplyr::mutate(DP = as.numeric(DP)) %>%
  ggplot(mapping = aes(x = DP, fill = SAMPLE)) +
  geom_histogram( alpha = 1, position = 'identity', binwidth = 50) +
  facet_grid(~SAMPLE, scales = 'free_x', space = 'free_x') +
  labs(title = 'distribution of total depth',
       x = 'depth of sequencing', y = 'number of variants') +
  theme_classic() + 
  my_theme +
  theme(strip.background = element_rect(fill = 'cornsilk2'),
        strip.text.x = element_text(size = 12))

#Q2 = plot_depth_dist + patchwork::plot_annotation(title = 'Q2')

#reads supporting the variant
plot_variant_depth_dist = final_annotation_table %>%
  dplyr::select(paste0('AO_', normal_sample_name), paste0('AO_', tumor_sample_name)) %>%
  tidyr::pivot_longer(cols = c(paste0('AO_', normal_sample_name), paste0('AO_', tumor_sample_name)),
                      names_to = 'SAMPLE', values_to = 'AO') %>%
  dplyr::mutate(AO = as.numeric(AO)) %>%
  ggplot(mapping = aes(x = AO, fill = SAMPLE)) +
  geom_histogram(alpha = 1, position = 'identity', binwidth = 50) +
  facet_grid(~SAMPLE, scales = 'free_x', space = 'free_x') +
  labs(title = 'distribution of total variant depth',
       x = 'depth of sequencing', y = 'number of variants') +
  theme_classic() + 
  my_theme +
  theme(strip.background = element_rect(fill = 'cornsilk2'),
        strip.text.x = element_text(size = 12))

#Q3 = plot_variant_depth_dist + patchwork::plot_annotation(title = 'Q3')

## Percentage of reads supporting the variant versus those supporting reference reads
plot_percAO_dist = final_annotation_table %>%
  dplyr::select(paste0('PERC_AO_', normal_sample_name), paste0('PERC_AO_', tumor_sample_name)) %>%
  tidyr::pivot_longer(cols = c(paste0('PERC_AO_', normal_sample_name), paste0('PERC_AO_', tumor_sample_name)),
                      names_to = 'SAMPLE', values_to = 'PERC_AO') %>%
  ggplot(mapping = aes(x = PERC_AO, fill = SAMPLE)) +
  geom_histogram() +
  facet_grid(~SAMPLE, scales = 'free_x', space = 'free_x') +
  labs(title = 'distribution of % reads supporting \nvariant and reference reads',
       x = '% reads', y = 'number of variants') +
  theme_classic() + 
  my_theme +
  theme(strip.background = element_rect(fill = 'cornsilk2'),
        strip.text.x = element_text(size = 12))

plot_percAO_normal0dist = final_annotation_table %>%
  dplyr::select(paste0('PERC_AO_', normal_sample_name), paste0('PERC_AO_', tumor_sample_name)) %>%
  dplyr::filter(!!as.symbol(paste0('PERC_AO_', normal_sample_name)) == 0) %>%
  tidyr::pivot_longer(cols = c(paste0('PERC_AO_', normal_sample_name), paste0('PERC_AO_', tumor_sample_name)),
                      names_to = 'SAMPLE', values_to = 'PERC_AO') %>%
  ggplot(mapping = aes(x = PERC_AO, fill = SAMPLE)) +
  geom_histogram() +
  facet_grid(~SAMPLE, scales = 'free_x', space = 'free_x') +
  labs(title = 'distribution of % reads supporting variant and \nreference reads when normal VAF = 0',
       x = '% reads', y = 'number of variants') +
  theme_classic() + 
  scale_x_continuous(limits = c(-2,10)) +
  my_theme +
  theme(strip.background = element_rect(fill = 'cornsilk2'),
        strip.text.x = element_text(size = 12))

#Q4 = (plot_percAO_dist / plot_percAO_normal0dist) + patchwork::plot_annotation(title = 'Q4')

# Allele frequency of variant from ExAC
plot_VAF_normalvaf0 = final_annotation_table %>%
  dplyr::select(paste0('FREQ_AO_', normal_sample_name), paste0('FREQ_AO_', tumor_sample_name), ExAC_ALL) %>%
  dplyr::mutate(ExAC_ALL = as.numeric(as.character(ExAC_ALL))) %>%
  dplyr::filter(!!as.symbol(paste0('FREQ_AO_', normal_sample_name)) == 0,
                !is.na(ExAC_ALL)) %>%
  dplyr::mutate(variant_id = row_number()) %>%
  tidyr::pivot_longer(cols = c(paste0('FREQ_AO_', normal_sample_name), paste0('FREQ_AO_', tumor_sample_name), 'ExAC_ALL'),
                      names_to = 'SAMPLE', values_to = 'FREQ_AO') %>%
  dplyr::mutate(FREQ_AO = as.numeric(FREQ_AO)) %>%
  ggplot(aes(x = SAMPLE, y = FREQ_AO)) + 
  geom_line(aes(group = factor(variant_id), color = factor(variant_id)), size = 1, show.legend = FALSE) + 
  geom_point(size = 0.1, color = 'red') +
  theme_classic() + 
  labs(title = 'AF for every mutation having normal VAF = 0',
       x = 'samples', y = 'AF') + 
  my_theme

#Q5 = plot_VAF_normalvaf0 + patchwork::plot_annotation(title = 'Q5')

final_plot = ((plot_variants_summary_box | plot_type_col) / (plot_variants_summary_col | plot_funcrefgene_col)) /
  (plot_depth_dist | plot_variant_depth_dist) /
  (plot_percAO_dist | plot_percAO_normal0dist) /
  plot_VAF_normalvaf0

ggsave(filename = paste0(file_path, '/', file_name, '.pdf'),
       plot = final_plot, device = 'pdf', 
       path = getwd(), scale = 1, width = 15, height = 20, 
       units = 'in', limitsize = TRUE, dpi = 100)
