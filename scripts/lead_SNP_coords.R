library(tidyverse)


# Find the coordinates of all lead SNPs in the library
# Output them in a table with SNP, chr, pos in hg37

# step 1 - make list of all lead SNPs in the library (filter out bladder and cervical)
# step 2 - are they all included in all-causal-SNP-16cancers.txt?
# step 3 - if not, output bed file and put through UCSC genome browser

res_merge <- read_tsv(here("output/res_merge_withhitannot_20230627.tsv"))
lead_dis <- res_merge %>%
  dplyr::select(lead_snp) %>%
  mutate(lead_snp = str_split(lead_snp, ",|;") ) %>%
  unnest(lead_snp)

lead_snps <- unique(lead_dis$lead_snp)
lead_snps <- lead_snps[!is.na(lead_snps)]

# Write out the SNPs
write(lead_snps, here("output/lead_snp_list.txt"))
# Go to https://genome.ucsc.edu/cgi-bin/hgTables
# genome: Human
# assembly: Feb 2009 (GRCh37/hg19)
# group: Variation
# track: dbSNP 155
# table: Common dbSNP (155) (dbSnp155Common)
# Upload list of lead SNPs for identifiers
# output format: BED
# output filename: data/lead_snp_coordinates_hg19_20230614






