## All R Code used in manuscript production
#load packages
library(dplyr)
library(ggplot2)
library(permute)
library(lattice)
library(vegan)
library(randomizeR)
library(ggpubr)
library(vctrs)
library(tidyverse)
library(reshape2)
library(rlang)
library(ggrepel)
library(plyr)
library(Maaslin2)
library(ggsci)
library(scales)
library(eulerr)
library(tidyr)
library(Hmisc)
library(phyloseq)
library(tidyjson)
library(readr)
library(grid)
library(stringr)
library(RColorBrewer)

## Randomization and distribution of casaes and controls across samples

# Load human milk sample aliquot ID information into R:
aliquots <- read.csv(file = 'V2 Final selected aliquots C1006 20220509.csv')
head(aliquots) # these are the samples I have for this project

# aliquots lists subjectnumber, group, PPID, Visit site, specimin label, ID sample

# subset samples into cases group and controls group
cases <- subset(aliquots, group=="Case")
controls <- subset(aliquots, group=="Control")

# add a column which is just numbering
cases$row_num <- seq.int(nrow(cases))
cases

controls$row_num <- seq.int(nrow(controls))
controls

# set seed to make randomization reproducible. 
set.seed(4256)

# randomize or shuffle the order of cases and controls.
data_cases <- cases[sample(1:nrow(cases)), ]
data_cases

data_controls <- controls[sample(1:nrow(controls)), ]
data_controls

# now I have randomly ordered cases and controls - These are called data_cases and data_controls

# Make a master list with these 

write.table(data_cases, "clipboard", sep="\t", row.names=FALSE)

write.table(data_controls, "clipboard", sep="\t", row.names=FALSE)

# copied both to excel. added a column of letters which ordered A, B, C... AA, BB etc. 1 per row for cases and 2 rows each letter for controls to randomize

cases_new <- read.csv(file = 'cases_random_plusletter.csv')
controls_new <- read.csv(file = 'controls_random_plusletter.csv')

total <- merge(cases_new, controls_new, by = "row_letter")

total <- rbind(cases_new, controls_new)
total2 <- arrange(total, by_group = row_letter)
total2

# now I have sorted the samples based on row letter. The samples are in an order case, control, control. Now I will add a numeric value to each sample.

total2$row_num <- seq.int(nrow(total2))
head(total2)

total2 <- subset(total2, select = c(subjectnumber, group, row_num))
head(total2)
colnames(total2) <- c('CHILD_ID', 'group', 'row_num')
# total2 lists in order (case, control, control) as they appear on the ELISA plates. labeled 1-300
# total2 sample order is used for all ELISAs so cases and controls are equally represented across plates


# Bring in total IgA ELISA results
ELISA_total <- read.csv(file = 'total_iga_results.csv')

# use this for matching subject number to CHILD ID
subjectnumber_to_ID <- read.csv(file = 'id_matching.csv')
head(subjectnumber_to_ID)
total <- merge(total2, ELISA_total, by = "CHILD_ID")
# now extract the total IgA results, CHILD_ID, and case/control group
total_iga <- subset(total, select = c(CHILD_ID, group, total_iga))

# Bring in other participant metadata from CHILD study
metadata <- read.csv('Metadata_BreastmilkIgA.csv')
metadata <- subset(metadata, select = c('SampleID', 'SubjectNumber', 'Visit', 'Blongum_clade', 'B.infantis.bin'))
colnames(metadata) <- c('SampleID', 'CHILD_ID', 'Visit', 'Blongum_clade', 'B.infantis.bin')
all_data <- metadata
shotgun_metadata <- read.delim('CP1006_filtered_simple_metadata_shotgun.txt')
head(shotgun_metadata)
shotgun_metadata <- subset(shotgun_metadata, select = c('SampleID', 'Exact_age_months'))
all_data <- merge(all_data, shotgun_metadata, by = 'SampleID')
# subset for 3mo data (the whole dataset includes 1y data too). 
# Create a master table with all data: rela_3mo
rela_3mo <- subset(all_data, Visit!='1 year')

# adding total IgA data
rela_3mo <- merge(rela_3mo, total_iga, by = 'CHILD_ID')


# adding all bacteria-specific IgA ELISA data

bbrev_iga <- read.csv('bbrev_iga_results.csv')
rela_3mo <- merge(rela_3mo, bbrev_iga, by = 'CHILD_ID')

para_iga <- read.csv('para_iga_results.csv')
rela_3mo <- merge(rela_3mo, para_iga, by = 'CHILD_ID')

blautia_iga <- read.csv("blautia_iga_results.csv")
rela_3mo <- merge(rela_3mo, blautia_iga, by = 'CHILD_ID')

dorea_iga <- read.csv('dorea_iga_results.csv')
colnames(dorea_iga) <- c('CHILD_ID', 'dorea_iga')
rela_3mo <- merge(rela_3mo, dorea_iga, by='CHILD_ID')

binf_iga <- read.csv('binf_iga_results.csv')
rela_3mo <- merge(rela_3mo, binf_iga, by = 'CHILD_ID')

strep_iga <- read.csv('strep_iga_results.csv')
rela_3mo <- merge(rela_3mo, strep_iga, by = 'CHILD_ID')

hhath_iga <- read.csv('child_hhath.csv')
rela_3mo <- merge(rela_3mo, hhath_iga, by = 'CHILD_ID')

rgnavus_iga <- read.csv('rgnavus_iga_results.csv')
rela_3mo <- merge(rela_3mo, rgnavus_iga, by = 'CHILD_ID')

bfrag_iga <- read.csv('bfrag_iga_results.csv')
rela_3mo <- merge(rela_3mo, bfrag_iga, by = 'CHILD_ID')

kpn_iga <- read.csv('kpne_results.csv')
rela_3mo <- merge(rela_3mo, kpn_iga, by = 'CHILD_ID')

entero_iga <- read.csv('enterococcus_results.csv')
rela_3mo <- merge(rela_3mo, entero_iga, by = 'CHILD_ID')

bdorei_iga <- read.csv('bdorei_iga_results.csv')
rela_3mo <- merge(rela_3mo, bdorei_iga, by = 'CHILD_ID')

eramosum_iga <- read.csv('eramosum_iga_results.csv')
rela_3mo <- merge(rela_3mo, eramosum_iga, by = 'CHILD_ID')

sparasan_iga <- read.csv('sparasanguinis_iga_results.csv')
rela_3mo <- merge(rela_3mo, sparasan_iga, by = 'CHILD_ID')

collinsella_iga <- read.csv('collinsella_iga_results.csv')
rela_3mo <- merge(rela_3mo, collinsella_iga, by = 'CHILD_ID')

buniformis_iga <- read.csv('buniformis_iga_results.csv')
rela_3mo <- merge(rela_3mo, buniformis_iga, by = 'CHILD_ID')

bbifidum_iga <- read.csv('bbifidum_iga_results.csv')
rela_3mo <- merge(rela_3mo, bbifidum_iga, by = 'CHILD_ID')

ecoli_iga <- read.csv('ecoli-results.csv')
ecoli_iga <- subset(ecoli_iga, select = c('CHILD_ID', 'ecoli_iga'))
rela_3mo <- merge(rela_3mo, ecoli_iga, by = 'CHILD_ID')

s_salivarius_iga <- read.csv('ssalivarius_iga_results.csv')
rela_3mo <- merge(rela_3mo, s_salivarius_iga, by = 'CHILD_ID')

cperfing_iga <- read.csv('cperfing_iga.csv')
rela_3mo <- merge(rela_3mo, cperfing_iga, by = 'CHILD_ID')

kvariicola_iga <- read.csv('kvariicola_iga_results.csv')
rela_3mo <- merge(rela_3mo, kvariicola_iga, by = 'CHILD_ID')

veillo_iga <- read.csv('child_veillo_results.csv')
rela_3mo <- merge(rela_3mo, veillo_iga, by = 'CHILD_ID')
head(rela_3mo)

# Adding IgA data normalized to total IgA

rela_3mo_norm <- rela_3mo %>% mutate(binf_norm = binf_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(bbrev_norm = bbrev_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(para_norm = para_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(dorea_norm = dorea_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(blautia_norm = blautia_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(strep_norm = strep_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(rgnavus_norm = rumino_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(bfrag_norm = bfrag_iga/total_iga/1000)
rela_3mo_norm$hhath_iga <- as.numeric(rela_3mo_norm$hhath_iga, na.rm = TRUE)
rela_3mo_norm <- rela_3mo_norm %>% mutate(hhath_norm = hhath_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(kpn_norm = kpn_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(entero_norm = entero_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(bdorei_norm = bdorei_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(sparasan_norm = sparasan_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(eramosum_norm = eramosum_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(collin_norm = collin_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(buinform_norm = buniform_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(ecoli_norm = ecoli_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(ssalivarius_norm = ssalivarius_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(kvariicola_norm = kvari_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(cperfing_norm = cperfing_iga/total_iga/1000)
rela_3mo_norm <- rela_3mo_norm %>% mutate(veillo_norm = veillo_iga/total_iga/1000)
head(rela_3mo_norm)
## Add in maternal allergy data:

# mom allergy metadata categories in metadata are:
# CHILD-97323 ever had asthma?
# CHILD-97324 asthma diagnosed by doc? 
# CHILD-97357 hay fever in last year?
# CHILD-97358 skin allergy symptoms in last year?
# CHILD-97364 food allergy symptoms in last year?

# read metadata file
output_data <- read.csv('output_data.csv')

output_data_mom_allergy <- subset(output_data, select = c('Participant_ID','CHILD.97323', 'CHILD.97324', 
                                                          'CHILD.97357', 'CHILD.97358', 'CHILD.97364'))
output_data_mom_allergy <- merge(output_data_mom_allergy, subjectnumber_to_ID, by = 'Participant_ID')
output_data_mom_allergy <- subset(output_data_mom_allergy, select = c('CHILD.97323', 
                                                                      'CHILD.97357', 'CHILD.97358', 
                                                                      'CHILD.97364', 'CHILD_ID'))
colnames(output_data_mom_allergy) <- c('ever_asthma', 'current_hay_fever', 
                                       'current_skin', 'current_food', 'CHILD_ID')
output_data_mom_allergy <- output_data_mom_allergy %>% mutate(atl_1_allergy = ifelse(ever_asthma == 1 |
                                                                                       current_hay_fever == 1 |
                                                                                       current_skin == 1 |
                                                                                       current_food == 1, 
                                                                                     "Y", "N"))

rela_3mo_norm <- merge(rela_3mo_norm, output_data_mom_allergy, by = 'CHILD_ID')
head(rela_3mo_norm)

## Erysipelotrichaceae data at family level (16s) at 3mo and 1y
# Erysipelo relative abundance at 3mo vs 1y
# 3mo: CHILD-150874
# 1y: CHILD-150902

output_data_eramosum <- subset(output_data, select = c('Participant_ID', 'CHILD.150874', 'CHILD.150902'))
output_data_eramosum <- merge(output_data_eramosum, subjectnumber_to_ID, by = 'Participant_ID')
output_data_eramosum$ery_3mo <- output_data_eramosum$CHILD.150874
output_data_eramosum$ery_1y <- output_data_eramosum$CHILD.150902
output_data_eramosum <- subset(output_data_eramosum, select = c('CHILD_ID', 'ery_3mo', 'ery_1y'))
head(output_data_eramosum)
# Remove NA values
output_data_eramosum <- na.omit(output_data_eramosum)

# Ensure the second and third columns are numeric
output_data_eramosum$ery_3mo <- as.numeric(output_data_eramosum$ery_3mo)
output_data_eramosum$ery_1y <- as.numeric(output_data_eramosum$ery_1y)
output_data_eramosum$deltaEry <- output_data_eramosum$ery_1y - output_data_eramosum$ery_3mo
rela_3mo_norm_eramosum <- merge(rela_3mo_norm, output_data_eramosum, by = 'CHILD_ID')

## Now I have two dataframes to work with:
# rela_3mo_norm has all IgA variable data
# rela_3mo_norm_eramosum has Erysipelotrichaceae data

setwd("~/OneDrive - UBC/Documents/Maaslin_IgA_ANALYSES")

# import microbiome metaddata from CHILD - subject #/ID/vist/age data
shotgun_metadata <- read.delim('CP1006_filtered_simple_metadata_shotgun.txt')
head(shotgun_metadata)
# subset for just 3mo data
sg_metadata_3m <- subset(shotgun_metadata, Visit == '3 month')

# import CHILD_ID/subject# data of my 300 samples
id_matching <- read.csv('id_matching.csv')

# subset metadata to just see samples that I used (subset of the 3 mo data)
my_3m <- merge(sg_metadata_3m, id_matching, by = 'subjectnumber')

# change the format of SampleID to match other data frames
my_3m$SampleID <- paste('X', my_3m$SampleID, sep = '')

# just select variables I want
sg_my_metadata_3m <- subset(my_3m, select = c('SampleID', 'Exact_age_months', 'site', 'CHILD_ID'))
head(sg_my_metadata_3m)
# import microbiome species data
shotgun_species <- read.delim('species.table.metaphlan_readcounts.txt')
# 637 species by 2952 samples
head(shotgun_species)

### Next, I need to convert species relative abundance data to a usable form
# from species data (in RPK form) get %RA 
# now need to only keep samples with more than 1000000 reads:
shotgun_species
row.names(shotgun_species) <- shotgun_species$clade_taxid
# take away phylogenetic classification info
phylogeny_table <- subset(shotgun_species, select = c('kingdom', 'phylum', 'class', 'order', 
                                                      'family', 'genus', 'species', 'clade_taxid'))
shotgun_species_condense <- subset(shotgun_species, select=-c(kingdom, phylum, class, order, family, genus, species))

## Filter data
#Remove columns where there are less than 1000000 reads
shotgun_species_filtered <- shotgun_species_condense[1 + which(colSums(shotgun_species_condense[-1]) > 1000000)]

#Remove rows with bugs less than 0.01% abundance
keep <- rowSums(shotgun_species_filtered) >= sum(colSums(shotgun_species_filtered))*0.0001 #only keep features with >0.01% of total abundance
shotgun_species_condense1 <- as.data.frame(subset(shotgun_species_filtered, keep == TRUE))
# 161 species kept
head(shotgun_species_condense1)

# add back clade_taxid 
taxav <- row.names(shotgun_species_condense1) # this gets the row numbers
taxa <- shotgun_species_condense[rownames(shotgun_species_condense) %in% taxav, ]
taxa <- subset(taxa, select = clade_taxid)
shotgun_species_condense1 <- cbind(taxa, shotgun_species_condense1)
# now I have a list of the samples which taxa filtered - called shotgun_species_condense1
taxa <- merge(taxa, phylogeny_table, by = 'clade_taxid')

# change to relative abundance
shotgun_species_RA <- shotgun_species_condense1 %>%
  mutate_at(vars(-1), funs(./sum(.))) 

shotgun_species_RA$clade_taxid <- rownames(shotgun_species_RA)
shotgun_species_RA_check <- merge(shotgun_species_RA, phylogeny_table, by = 'clade_taxid')

# check work:
# column 2, entry 27 is 0.0796383
# column 2 total:
colSums(shotgun_species_condense[2])
# 2124836
# column 2 entry 26 is 157161
shotgun_species_condense[27,2]
157161/2124836
taxav <- row.names(shotgun_species_condense1)
shotgun_ra_filtered <- shotgun_species_RA[rownames(shotgun_species_RA) %in% taxav, ]
# now I have converted to RA and then filtered for top bugs and added species names back to data frame

# now need to select just my 3mo samples using metadata 
all_3m_metadata <- IgA_data
subjects <- all_3m_metadata$SampleID
shotgun_ra_3m <- subset(shotgun_ra_filtered, select = subjects)

# now I have only 3m samples, filtered first based on read count, then species were filtered based on abundance
# above 0.01%. then just selected my 288 samples for analysis

# Now change clade_taxid to species
species <- subset(phylogeny_table, select = c('species', 'clade_taxid'))
head(species)
shotgun_ra_3m$clade_taxid <- row.names(shotgun_ra_3m)
shotgun_ra_3m <- merge(species, shotgun_ra_3m, by = 'clade_taxid')
shotgun_ra_3m_species <- subset(shotgun_ra_3m, select=-clade_taxid)
head(shotgun_ra_3m_species)
my_3m_species <- shotgun_ra_3m_species

## Now I want to determine which species are most abundant and prevalent in my sample of infants:
head(phylogeny_table)
phylogeny_table <- phylogeny_table[-8]
row.names(phylogeny_table) <- phylogeny_table$species
tax<-tax_table(as.matrix(phylogeny_table))
tax
head(my_3m_species)
row.names(my_3m_species) <- my_3m_species$species
my_3m_species <- my_3m_species[-1]
otu<-otu_table(my_3m_species, taxa_are_rows = TRUE)

row.names(all_3m_metadata) <- all_3m_metadata$SampleID
sample<-sample_data(all_3m_metadata)
sample_names(otu)
sample_names(tax)
sample_names(sample)
#merge into phyloseq
physeq.3m = phyloseq(otu,tax,sample)
physeq.3m

## From this phyloseq object, I want to get top 10 most abundant families:
# Check the aggregated family-level object
family<-tax_glom(physeq.3m, taxrank = "family", NArm=FALSE)
family.3m<-family@otu_table@.Data
rownames(family.3m)<-phylogeny_table$family[match(rownames(family.3m), rownames(phylogeny_table))]
family

otu_table(family)
tax_table(family)
sort(taxa_sums(family))

fam_3m <- prune_taxa(
  names( sort( taxa_sums(family), decreasing = T )[ 1:30] ),
  family)
fam_3m
sample_data(fam_3m)$SampleType <- factor(sample_names(fam_3m))


fam_3m_df<- psmelt(fam_3m)
head(fam_3m_df)

nb.cols <- 25
mycolors<- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
fam_3m_df %>% ggplot(aes(y = Abundance, x = reorder(family, -Abundance))) + 
  scale_fill_manual(values = mycolors) +
  geom_boxplot(aes(fill = fam_3m_df$family), show.legend = F)+
  #facet_wrap(~ fam_3m_df$group, scales = "free_x") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.1)) + 
  labs(x = "Family", y = "Relative Abundance")
# This graph shows the top 10 most abundant families in the infants 
# families:               
# "Bifidobacteriaceae" 
# "Enterobacteriaceae" 
# "Bacteroidaceae"     
# "Lachnospiraceae"    
# "Enterococcaceae"    
# "Erysipelotrichaceae"
# "Clostridiaceae"     
# "Streptococcaceae"   
# "Tannerellaceae"     
# "Coriobacteriaceae"  
### NEED TO FIGURE OUT THIS

## within these, I want to see the species with prevalence > 50
my_3m_species_p <- my_3m_species
head(my_3m_species_p)
my_3m_species_p <- mutate_all(my_3m_species_p, function(x) as.numeric(as.character(x)))
my_3m_species_prev <- my_3m_species_p %>% mutate_all(funs(ifelse(.>0,1,0)))
my_3m_species_prev$sum <- rowSums(my_3m_species_prev)
my_3m_species_prev$species <- row.names(my_3m_species_prev)
my_3m_species_prev <- merge(my_3m_species_prev, phylogeny_table, by = 'species')
head(my_3m_species_prev)
row.names(my_3m_species_prev) <- my_3m_species_prev$species
# now I have a data frae with 1 or 0 for infants who have or do not have each species
# now for each of top 10 families, find most prevalent bugs:
# Bifidobacteriaceae
Bifido <- subset(my_3m_species_prev, family == 'Bifidobacteriaceae')
row.names(Bifido) <- Bifido$species
Bifido <- subset(Bifido, select=-c(species, kingdom, phylum, class, order, family, genus))
Bifido <- subset(Bifido, select = sum)
Bifido
# top are B longum, B breve, B bifidum, B pseudocatenulatum
# prev are B. longum 162, B breve 94, B bifidum 66, B dentium 56, pseudocatenulatum 43, animalis 31

# Enterobacteriaceae
Entero <- subset(my_3m_species_prev, family == 'Enterobacteriaceae')
row.names(Entero) <- Entero$species
Entero <- subset(Entero, select=-c(species, kingdom, phylum, class, order, family, genus))
Entero <- subset(Entero, select = sum)
Entero
# top are E coli, K pneumonieae, K oxytoca, K variicola
# prev are E. coli 201, K pneumoniae 77,K variicola 66, K quasipneumoniae 59, K oxycota 50, K michigansis 40

# Bacteroidaceae
Bacteroid <- subset(my_3m_species_prev, family == 'Bacteroidaceae')
row.names(Bacteroid) <- Bacteroid$species
Bacteroid <- subset(Bacteroid, select=-c(species, kingdom, phylum, class, order, family, genus))
Bacteroid <- subset(Bacteroid, select = sum)
Bacteroid
# top are B fragilis, B uniformis, B vulgatus, B dorei, B theta
# prev are B fragilis 57, B uniformis 56, B vulgatus 49, B ovatus 35, theta 30

# Lachnospiracaeae
Lachno <- subset(my_3m_species_prev, family == 'Lachnospiraceae')
row.names(Lachno) <- Lachno$species
Lachno <- subset(Lachno, select=-c(species, kingdom, phylum, class, order, family, genus))
Lachno <- subset(Lachno, select = sum)
Lachno
# R. gnavus, B. wexlerae, F. saccharivorans, Roseburia intestinalis, Roseburia faecis
# prev are R. gnavus 73, B wex 33, C symbiosum 18

# Enterococcaceae
Enterococ <- subset(my_3m_species_prev, family == 'Enterococcaceae')
row.names(Enterococ) <- Enterococ$species
Enterococ <- subset(Enterococ, select=-c(species, kingdom, phylum, class, order, family, genus))
Enterococ <- subset(Enterococ, select = sum)
Enterococ
# Enterococcus faecalis, Enterococcus durans, avium, casseliflavus
# prev are Enterococcus faecalis 181, Enterococcus avium 22

# Erysipelotrichaceae
Erysipo <- subset(my_3m_species_prev, family == 'Erysipelotrichaceae')
row.names(Erysipo) <- Erysipo$species
Erysipo <- subset(Erysipo, select=-c(species, kingdom, phylum, class, order, family, genus))
Erysipo <- subset(Erysipo, select = sum)
Erysipo
# Erysipelatoclostridium ramosum, C innocuum
# top are E ramosum 80, C innocum 32

# Clostridicaceae
Clostrid <- subset(my_3m_species_prev, family == 'Clostridiaceae')
row.names(Clostrid) <- Clostrid$species
Clostrid <- subset(Clostrid, select=-c(species, kingdom, phylum, class, order, family, genus))
Clostrid
Clostrid <- subset(Clostrid, select = sum)
Clostrid
# C. neonatale, H hathewayi, C perfingens
# prev are C perfingens 45, C neonatale 40, C paraputrificum 28, HH 29

# Streptococcacaeae
Strep <- subset(my_3m_species_prev, family == 'Streptococcaceae')
row.names(Strep) <- Strep$species
Strep <- subset(Strep, select=-c(species, kingdom, phylum, class, order, family, genus))
Strep <- subset(Strep, select = sum)
Strep
# Strep parasanguinis, Strep salivarius, Strep vestibularis, Strep lutetiensis, no S gallo...
# prev are S salivarius 148, S parasangunis 146, S mitis 121, S vestibulus 39, S thermo 33

# Tannerellaceae
Tanner <- subset(my_3m_species_prev, family == 'Tannerellaceae')
row.names(Tanner) <- Tanner$species
Tanner <- subset(Tanner, select=-c(species, kingdom, phylum, class, order, family, genus))
Tanner <- subset(Tanner, select = sum)
Tanner
# Parabacteroids distasonis was highest (prev 59)

# Coriobacteraceae
Cori <- subset(my_3m_species_prev, family == 'Coriobacteriaceae')
row.names(Cori) <- Cori$species
Cori <- subset(Cori, select=-c(species, kingdom, phylum, class, order, family, genus))
Cori <- subset(Cori, select = sum)
Cori
# none with prev>50

# From this list, I chose all species with prevalence > 50, or for Strep just top 2 species
# THESE ARE LISTED IN TABLE 1
# These were the bacteria grown and tested in IgA-bacteria ELISAs

# To determine average abundance of these among infants who have them:
my_3m_species_a <- my_3m_species
head(my_3m_species_a)
my_3m_species_a <- mutate_all(my_3m_species_a, function(x) as.numeric(as.character(x)))

my_3m_species_a$average <- apply(my_3m_species_a, 1, function(x) {
  non_zero <- as.numeric(x[x != 0])
  if (length(non_zero) == 0) {
    return(0)  # or NA if preferred
  } else {
    return(mean(non_zero))
  }
})

# as well as range in abundance:
my_3m_species_a$range <- apply(my_3m_species_a, 1, function(x) {
  non_zero <- as.numeric(x[x != 0])
  if (length(non_zero) == 0) {
    return(0)  # or NA if you prefer
  } else {
    return(max(non_zero) - min(non_zero))
  }
})


my_3m_species_a$species <- row.names(my_3m_species_a)
my_3m_species_ab_range <- subset(my_3m_species_a, select = c('species', 'average', 'range'))
head(my_3m_species_ab_range)

# average for each species on list:
my_3m_species_a <- subset(my_3m_species_a, select = c('species', 'average'))
head(my_3m_species_a)


# now I have a data frae with 1 or 0 for infants who have or do not have each species
# now for each of top 10 families, find most prevalent bugs:
# Bifidobacteriaceae
Bifido <- subset(my_3m_species_prev, family == 'Bifidobacteriaceae')
row.names(Bifido) <- Bifido$species
Bifido <- subset(Bifido, select=-c(species, kingdom, phylum, class, order, family, genus))
Bifido <- subset(Bifido, select = sum)
Bifido

## First I want to look at general patterns in IgA data:
head(rela_3mo_norm)
IgA_data <- rela_3mo_norm
head(IgA_data)
iga_data_justvariables <- subset(IgA_data, select = c(bbrev_iga, para_iga, bfrag_iga, binf_iga, 
                                                      rumino_iga, kpn_iga, entero_iga,  eramosum_iga, 
                                                      sparasan_iga, buniform_iga, bbifidum_iga, ecoli_iga, ssalivarius_iga, kvari_iga, cperfing_iga, veillo_iga))

colnames(iga_data_justvariables) <- c('Bifidobacterium breve','Parabacteroides distasonis', 'Bacteroides fragilis', 
                                      'Bifidobacterium longum infantis',
                                      'Ruminococcus gnavus', 'Klebsiella pneumoniae', 
                                      'Enterococcus faecalis', 'Erysipelatoclostridium ramosum',
                                      'Streptococcus parasanguinis', 'Bacteroides uniformis', 
                                      'Bifidobacterium bifidum', 'Escherichia coli', 'Streptococcus salivarius', 'Klebsiella variicola',
                                      'Clostridium perfingens')
numeric_data_long <- gather(iga_data_justvariables, key = "Variable", value = "Value")
head(numeric_data_long)

# preparing a graph of just the variation in specific IgA across bugs (Figure 1B)
numeric_data_long = na.omit(numeric_data_long)
numeric_data_long <- subset(numeric_data_long, Value < 60000) # getting rid of outliers
numeric_data_long$Value <- numeric_data_long$Value/1000
numeric_data_long$Variable = reorder(numeric_data_long$Variable, numeric_data_long$Value, mean)
iga_levels <- ggplot(numeric_data_long, aes(x = Variable, y = Value)) +
  geom_boxplot(fill = 'lightskyblue2') +
  theme_pubclean() +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2, color = "maroon") +
  labs(title = "Average quantity of binding milk SIgA",
       x = "tested species",
       y = "milk IgA binding capacity (ug/mL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 'italic'))
iga_levels
ggsave('iga-binding.svg', plot = iga_levels)
# THIS IS FIGURE 1B

# Now I need to get IgA data into a format that can be used for comparison to species RA data
head(IgA_data)
# to make sample IDs match:
IgA_data$SampleID <- paste("X", IgA_data$SampleID, sep = '')
IgA_data <- subset(IgA_data, select=-CHILD_ID)
all_3m_metadata  <- IgA_data
row.names(all_3m_metadata) <- all_3m_metadata$SampleID

## species analysis - run Maaslin2 on IgA variables + species relative abundance
# this will look for associations between bugs and each variable
# input data: the species RA data frame
# input metadata: all my metadata effects

# exact age (months) of infant:
head(all_3m_metadata)
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_exact_age",
  max_significance = 0.6,
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = "Exact_age_months")
exact_age_results <- read_tsv('output_exact_age/all_results.tsv')
exact_age_results$feature <- gsub("_", ' ', exact_age_results$feature)


# Specify plot_maaslin function:
plot_maaslin <- function(df) {
  # Plots the average, min, and max inflammation over time.
  # Input is character string of a csv file.
  features.keep <- subset(df)$feature
  df.2 <- subset(df,feature %in% features.keep)
  df.2$sig <- ifelse(df.2$pval<0.05,1,0) 
  ggplot(df.2, aes(x = coef, y = reorder(feature, coef), color=value, alpha=sig)) +
    geom_vline(aes(xintercept = 0), size = 1, linetype = "dashed") + #add line at 1 odds ratio
    geom_point(size=2, position = position_dodge(width = 0.75)) +
    geom_errorbarh(aes(xmin = coef + stderr, xmax = coef - stderr),size=1,
                   position = position_dodge(width = 0.75))+
    theme_bw()+
    #geom_text(aes(label = ifelse(sig==1,paste0("FDR=", round(qval,3)),"")), vjust=-1)+
    #coord_cartesian(xlim = c(0.01,3))+
    ylab("") +
    guides(alpha=FALSE)+
    xlab("Coefficient")+
    #ggtitle("Case vs Control")+
    theme(
      legend.position = "none",
      axis.text.y=element_text(face="italic",size=10), 
      #axis.title.y=element_text(face="bold",size=14), 
      #axis.title.x =element_blank(),
      plot.title = element_text(face="bold", size=16), 
      strip.text = element_text(face="bold",size=14), 
      panel.border = element_rect(fill=NA, colour = "black", size=2),
      panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "transparent",colour = NA), 
      plot.background = element_rect(fill = "transparent",colour = NA)) }

exact_age_vs_bacteria <- plot_maaslin(exact_age_results)
ggsave('exact_age_vs_bacteria.svg', plot = exact_age_vs_bacteria)
# THIS IS FIGURE 2A




# total IgA:
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_total_iga",
  max_significance = 0.6,
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("total_iga", "Exact_age_months"))
total_iga_results <- read_tsv('output_total_iga/all_results.tsv')
total_iga_results <- subset(total_iga_results, metadata == 'total_iga')

# Species specific IgA ELISA data associations
# B. longum
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_binf_iga",
  normalization = "NONE",
  transform = "LOG",
  max_significance = 0.6,
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("binf_iga", "Exact_age_months"))

# B. breve
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_bbrev_iga",
  normalization = "NONE",
  transform = "LOG",
  max_significance = 0.6,
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("bbrev_iga", "Exact_age_months"))

# B. bifidum
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_bbif_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("bbifidum_iga", "Exact_age_months"))

# E. coli
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_ecoli_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("ecoli_iga", "Exact_age_months"))

# K. pneumoniae
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_kpn_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("kpn_iga", "Exact_age_months"))

# K. variicola
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_kvari_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("kvari_iga", "Exact_age_months"))

# B. fragilis
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_bfrag_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("bfrag_iga", "Exact_age_months"))

# B. uniformis
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_buniform_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("buniform_iga", "Exact_age_months"))

# R. gnavus
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_rumino_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("rumino_iga", "Exact_age_months"))

# E. faecalis
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_entero_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("entero_iga", "Exact_age_months"))

# E. ramosum
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_eramosum_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("eramosum_iga", "Exact_age_months"))

# C. perfingens
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_cperfing_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("cperfing_iga", "Exact_age_months"))

# S. salivarius
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_ssalivarius_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("ssalivarius_iga", "Exact_age_months"))

# S. parasanguinis
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_sparasan_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("sparasan_iga", "Exact_age_months"))

# P. distasonis
fit_data = Maaslin2(
  input_data = my_3m_species, 
  input_metadata = all_3m_metadata,
  output = "output_para_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  max_significance = 0.5,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("para_iga", "Exact_age_months"))

## START HERE
# need to do them all individually

binf_iga_results <- read_tsv('output_binf_iga/all_results.tsv')
binf_iga_results <- subset(binf_iga_results, metadata == 'binf_iga')
bbrev_iga_results <- read_tsv('output_bbrev_iga/all_results.tsv')
bbrev_iga_results<- subset(bbrev_iga_results, metadata == 'bbrev_iga')
bbif_iga_results <- read_tsv('output_bbif_iga/all_results.tsv')
bbif_iga_results <- subset(bbif_iga_results, metadata == 'bbifidum_iga')
ecoli_iga_results <- read_tsv('output_ecoli_iga/all_results.tsv')
ecoli_iga_results <- subset(ecoli_iga_results, metadata == 'ecoli_iga')
kpn_iga_results <- read_tsv('output_kpn_iga/all_results.tsv')
kpn_iga_results <- subset(kpn_iga_results, metadata == 'kpn_iga')
kvari_iga_results <- read_tsv('output_kvari_iga/all_results.tsv')
kvari_iga_results <- subset(kvari_iga_results, metadata == 'kvari_iga')
bfrag_iga_results<- read_tsv('output_bfrag_iga/all_results.tsv')
bfrag_iga_results<- subset(bfrag_iga_results, metadata == 'bfrag_iga')
buniform_iga_results <- read_tsv('output_buniform_iga/all_results.tsv')
buniform_iga_results <- subset(buniform_iga_results, metadata == 'buniform_iga')
rumino_iga_results<- read_tsv('output_rumino_iga/all_results.tsv')
rumino_iga_results<- subset(rumino_iga_results, metadata == 'rumino_iga')
entero_iga_results<- read_tsv('output_entero_iga/all_results.tsv')
entero_iga_results<- subset(entero_iga_results, metadata == 'entero_iga')
eramosum_iga_results<- read_tsv('output_eramosum_iga/all_results.tsv')
eramosum_iga_results<- subset(eramosum_iga_results, metadata == 'eramosum_iga')
cperfing_iga_results<- read_tsv('output_cperfing_iga/all_results.tsv')
cperfing_iga_results<- subset(cperfing_iga_results, metadata == 'cperfing_iga')
ssalivarius_iga_results<- read_tsv('output_ssalivarius_iga/all_results.tsv')
ssalivarius_iga_results<- subset(ssalivarius_iga_results, metadata == 'ssalivarius_iga')
sparasan_iga_results<- read_tsv('output_sparasan_iga/all_results.tsv')
sparasan_iga_results<- subset(sparasan_iga_results, metadata == 'sparasan_iga')
para_iga_results<- read_tsv('output_para_iga/all_results.tsv')
para_iga_results<- subset(para_iga_results, metadata == 'para_iga')

# individuate each data frame so that all can be merged
species_binf <- subset(binf_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_binf) <- c('feature', 'binf_coef', 'binf_pval', 'binf_stderr')
species_bbif <- subset(bbif_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_bbif) <- c('feature', 'bbif_coef', 'bbif_pval', 'bbif_stderr')
species_bbrev <- subset(bbrev_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_bbrev) <- c('feature', 'bbrev_coef', 'bbrev_pval', 'bbrev_stderr')
species_ecoli <- subset(ecoli_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_ecoli) <- c('feature', 'ecoli_coef', 'ecoli_pval', 'ecoli_stderr')
species_kpn <- subset(kpn_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_kpn) <- c('feature', 'kpn_coef', 'kpn_pval', 'kpn_stderr')
species_kvari <- subset(kvari_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_kvari) <- c('feature', 'kvari_coef', 'kvari_pval', 'kvari_stderr')
species_bfrag <- subset(bfrag_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_bfrag) <- c('feature', 'bfrag_coef', 'bfrag_pval', 'bfrag_stderr')
species_buniform <- subset(buniform_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_buniform) <- c('feature', 'buniform_coef', 'buniform_pval', 'buniform_stderr')
species_rumino <- subset(rumino_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_rumino) <- c('feature', 'rumino_coef', 'rumino_pval', 'rumino_stderr')
species_entero <- subset(entero_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_entero) <- c('feature', 'entero_coef', 'entero_pval', 'entero_stderr')
species_eramosum <- subset(eramosum_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_eramosum) <- c('feature', 'eramosum_coef', 'eramosum_pval', 'eramosum_stderr')
species_cperfing <- subset(cperfing_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_cperfing) <- c('feature', 'cperfing_coef', 'cperfing_pval', 'cperfing_stderr')
species_ssalivarius <- subset(ssalivarius_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_ssalivarius) <- c('feature', 'ssalivarius_coef', 'ssalivarius_pval', 'ssalivarius_stderr')
species_sparasan <- subset(sparasan_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_sparasan) <- c('feature', 'sparasan_coef', 'sparasan_pval', 'sparasan_stderr')
species_para <- subset(para_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_para) <- c('feature', 'para_coef', 'para_pval', 'para_stderr')
species_total <- subset(total_iga_results, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(species_total) <- c('feature', 'total_coef', 'total_pval', 'total_stderr')

# put all together in a table:

species_table <- merge(species_bfrag, species_bbrev, by = 'feature')
head(species_table)
species_table <- merge(species_table, species_binf, by = 'feature')
species_table <- merge(species_table, species_bbif, by = 'feature')
species_table <- merge(species_table, species_buniform, by = 'feature')
species_table <- merge(species_table, species_ecoli, by = 'feature')
species_table <- merge(species_table, species_kpn, by = 'feature')
species_table <- merge(species_table, species_kvari, by = 'feature')
species_table <- merge(species_table, species_sparasan, by = 'feature')
species_table <- merge(species_table, species_ssalivarius, by = 'feature')
species_table <- merge(species_table, species_rumino, by = 'feature')
species_table <- merge(species_table, species_entero, by = 'feature')
species_table <- merge(species_table, species_eramosum, by = 'feature')
species_table <- merge(species_table, species_cperfing, by = 'feature')
species_table <- merge(species_table, species_para, by = 'feature')
species_table <- merge(species_table, species_total, by = 'feature')

species_table_coef <- subset(species_table, select = c(feature, binf_coef, bbrev_coef, bbif_coef, 
                                                       bfrag_coef, buniform_coef, ecoli_coef, kpn_coef, 
                                                       kvari_coef, sparasan_coef, ssalivarius_coef, rumino_coef, 
                                                       entero_coef, eramosum_coef, cperfing_coef, para_coef, total_coef))
colnames(species_table_coef) <- c('feature',  'B. longum IgA', 'B. breve IgA', 'B. bifidum IgA', 
                                  'B. fragilis IgA', 'B. uniformis IgA', 'E. coli IgA', 'K. pneumoniaeae IgA', 
                                  'K. variicola IgA', 'S. parasanguinis IgA', 'S. salivarius IgA', 
                                  'R. gnavus IgA', 'E. faecalis IgA', 'E. ramosum IgA', 'C. perfingens IgA',
                                  'P. disasonis IgA', 'total IgA')
head(species_table_coef)
data_melt <- melt(species_table_coef)
head(data_melt)


# now want to get all species-specific IgA associations and total IgA associations for a bubble plot (Figure 1C)
colnames(species_bbrev) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_bbrev <- subset(species_bbrev, feature == 'Bifidobacterium_breve')
iga_bbrev <- as.data.frame(iga_bbrev)
colnames(species_bbif) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_bbifidum <- subset(species_bbif, feature == 'Bifidobacterium_bifidum')
iga_bbifidum <- as.data.frame(iga_bbifidum)
colnames(species_binf) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_binf <- subset(species_binf, feature == 'Bifidobacterium_longum')
iga_binf <- as.data.frame(iga_binf)
colnames(species_ecoli) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_ecoli <- subset(species_ecoli, feature == 'Escherichia_coli')
iga_ecoli <- as.data.frame(iga_ecoli)
colnames(species_kpn) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_kpn <- subset(species_kpn, feature == 'Klebsiella_pneumoniae')
iga_kpn <- as.data.frame(iga_kpn)
colnames(species_kvari) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_kvari <- subset(species_kvari, feature == 'Klebsiella_variicola')
iga_kvari <- as.data.frame(iga_kvari)
colnames(species_bfrag) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_bfrag <- subset(species_bfrag, feature == 'Bacteroides_fragilis')
iga_bfrag <- as.data.frame(iga_bfrag)
colnames(species_buniform) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_buniform <- subset(species_buniform, feature == 'Bacteroides_uniformis')
iga_buniform <- as.data.frame(iga_buniform)
colnames(species_rumino) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_rumino <- subset(species_rumino, feature == 'Ruminococcus_gnavus')
iga_rumino <- as.data.frame(iga_rumino)
colnames(species_entero) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_entero <- subset(species_entero, feature == 'Enterococcus_faecalis')
iga_entero <- as.data.frame(iga_entero)
colnames(species_eramosum) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
iga_eramosum <- subset(species_eramosum, feature == 'Erysipelatoclostridium_ramosum')
iga_eramosum <- as.data.frame(iga_eramosum)
head(species_cperfing)
colnames(species_cperfing) <- c('feature', 'spec_coef', 'spec_stderr')
stacked_df_con <- stacked_df_con %>% mutate(spec_sig = ifelse(abs(spec_pval) < 0.1, 'y', 'n'))
stacked_df_con <- stacked_df_con %>% mutate(total_sig = ifelse(abs(total_pval) < 0.1, 'y', 'n'))

stacked_df_con <- subset(stacked_df_con, select = c(feature, spec_coef, spec_size, spec_sig, total_coef, total_size, total_sig))
# Reshape the dataframe for 'spec'

df_spec <- stacked_df_con %>%
  select(feature, spec_coef, spec_size, spec_sig) %>%
  mutate(category = "spec", coef = spec_coef, size = spec_size, sig = spec_sig) %>%
  select(feature, category, coef, size, sig)

# Reshape the dataframe for 'total'
df_total <- stacked_df_con %>%
  select(feature, total_coef, total_size, total_sig) %>%
  mutate(category = "total", coef = total_coef, size = total_size, sig = total_sig) %>%
  select(feature, category, coef, size, sig)
=
  select(feature, category, coef, size_abs, sig)

df_combined
df_combined$category <- factor(df_combined$category, levels = c("total", "spec"))
df_combined <- df_combined %>%
  mutate(category = recode(category, "spec" = "specific IgA", "total" = "total IgA"))

# now want to create a bubble plot:
# add a column with an asterisk if significant
df_combined$asterisk <- ifelse(df_combined$sig == 'y', "*", "")
df_combined
df_combined$feature <- gsub("_", " ", df_combined$feature)
p <- ggplot(df_combined, aes(x=category, y=feature, size=abs(coef))) + 
  geom_point(shape=21, aes(fill=size_abs)) + 
  theme(axis.text.y = element_text(face = "italic")) +
  geom_rect(aes(xmin=2, xmax=Inf, ymin=-Inf, ymax=Inf), fill="gray90", alpha=0.5, size = 0.1) +
  geom_text(aes(label=asterisk), vjust=0.72, size=8, color = 'white') 
p
legend_labels <- c("y_neg" = "(-) coef > 0.3", "n_neg" = "(-)", "y_pos" = "(+) coef > 0.3", "n_pos" = "(+)")

p <- p + scale_fill_manual(values = c("y_neg" = "indianred3", "n_neg" = "indianred1", "y_pos" = "deepskyblue4", "n_pos" = "deepskyblue1"),
                           breaks = names(legend_labels), labels = legend_labels)
p <- p + scale_colour_manual(values=c("grey80", "black"))
p <- p + scale_size_continuous(range = c(0.5,15))
p <- p + theme(axis.text.x  = element_blank(), axis.title.y = element_blank(), legend.key.size=unit(0.75, "cm"))
p <- p + scale_y_discrete(limits = rev(levels(df_combined$feature)))
#p <-  p + guides(color = guide_legend(override.aes = list(size=6)))
p <- p + facet_grid(.~category, scales = "free", space = "free") +
  labs(fill = "Direction & Size", color = "Direction & Significance", size = "Absolute Coefficient") 

p <- p + theme(strip.background = element_rect(fill="gray85"),
               panel.background = element_rect(fill="white"),
               panel.border = element_rect(colour="black", linetype="solid", fill="transparent") 
)
dev.new(width=8, height=13)
p
ggsave('total_and_specific_iga_associations.svg', plot = p)
# this bubble plot is Figure 1C


# Associations between IgA binding variables, Figure S1A:
head(all_3m_metadata)
iga_data_justvariables <- subset(all_3m_metadata, select = c("total_iga", "bbrev_iga", "binf_iga", "bbifidum_iga",
                                                              "bfrag_iga", "buniform_iga", "ecoli_iga", "kpn_iga",
                                                              "kvari_iga", "rumino_iga", "entero_iga",
                                                              "eramosum_iga", "cperfing_iga", "ssalivarius_iga",
                                                              "sparasan_iga", "para_iga"))
colnames(iga_data_justvariables) <- c("Total", "B. breve", "B. longum infantis", "B. bifidum",
                                           "B, fragilis", "B. uniformis", "E. coli", "K. pneumonieae",
                                           "K. variicola", "R. gnavus", "E. faecalis",
                                           "E. ramosum", "C. perfringens", "S. salivarius",
                                          "S. parasanguinis", "P. distansonis")

# Calculate Pearson correlation matrix
cor_matrix <- cor(iga_data_justvariables, method = "pearson", use = "complete.obs")


# Perform hierarchical clustering on the rows and columns of the correlation matrix
row_order <- hclust(dist(cor_matrix))$order
col_order <- hclust(dist(t(cor_matrix)))$order

# Reorder the correlation matrix
cor_matrix_ordered <- cor_matrix[row_order, col_order]

# Melt the reordered correlation matrix for ggplot
melted_cor_matrix <- melt(cor_matrix_ordered)

# Custom function to make labels italics (same as before)
custom_labels <- function(labels) {
  sapply(labels, function(label) {
    if (label == "Total") label else bquote(italic(.(label)))
  })
}

# Plot the heatmap with reordered correlation matrix
correlation_iga_binding_variables <- ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab", name="Pearson\nCorrelation") +
  labs(x = NULL, y = NULL) +  # Remove axis labels
  scale_x_discrete(labels = custom_labels, limits = colnames(cor_matrix_ordered)) +  # Apply custom labels for x-axis
  scale_y_discrete(labels = custom_labels, limits = rownames(cor_matrix_ordered)) +  # Apply custom labels for y-axis
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed()

# Display the heatmap
print(correlation_iga_binding_variables)
ggsave('correlation_iga_binding_variables.svg', correlation_iga_binding_variables)
# THIS IS FIGURE S1A


##  Total SIgA by infant and maternal atopy status (Supplemental Figure 1)
infant_allergy_iga <- ggplot(all_3m_metadata, aes(x = group, y = total_iga)) +
  geom_boxplot(color = "Deep sky blue", size = 1.2, outlier.shape = NA) +
  geom_jitter(width = 0.35, size = 2, color = "black") +
  theme_pubclean() +
  labs(x = "Infant Group", y = "Total SIgA (g/L) in milk", 
       title = "Total SIgA by infant atopy") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("Case", "Control")),
                     label.y = 5) +  # Set label position above the y-axis limit 
  coord_cartesian(ylim = c(0, 6))  # Set y-axis limits without removing data points
infant_allergy_iga
ggsave('infant_allergy_iga.svg', infant_allergy_iga)
# THIS IS FIGURE S1B

# Getting mom atopy data from spt dataset from CHILD Study:
# adding in spt data:

spt_data <- read.csv('spt_data_mod.csv')
head(spt_data)
spt_data_mom <- subset(spt_data, select = c('SampleID', 'mom_atopy'))
head(spt_data_mom)
# this gives a list of mom atopy status (0 or 1) by Sample ID
spt_data_mom$SampleID <- paste("X", spt_data_mom$SampleID, sep = '')
head(spt_data_mom)
all_3m_metadata <- merge(all_3m_metadata, spt_data_mom, by = 'SampleID')
class(all_3m_metadata$mom_atopy)
all_3m_metadata$mom_atopy <- ifelse(all_3m_metadata$mom_atopy == 1, "atopic", "non-atopic")
mom_allergy_iga <- ggplot(all_3m_metadata, aes(x = mom_atopy, y = total_iga)) +
  geom_boxplot(color = "Deep sky blue", size = 1.2, outlier.shape = NA) +
  geom_jitter(width = 0.35, size = 2, color = "black") +
  theme_pubclean() +
  labs(x = "Maternal group", y = "Total SIgA (g/L) in milk", 
       title = "Total milk SIgA by maternal atopy status") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("non-atopic", 'atopic')),
                     label.y = 4) +  # Set label position above the y-axis limit 
  coord_cartesian(ylim = c(0, 6))  # Set y-axis limits without removing data points
mom_allergy_iga
ggsave('mom_allergy_iga.svg', mom_allergy_iga)
# THIS IS FIGURE S1C

## Olink Biomarker Data (Supplemental Figure 1)

# Import Olink Biomarker data
olink_data <- read.csv('olink-data.csv')
colnames(olink_data)[1] <- "SampleID"
# import metadata (1y and 5y and exat age)
olink_guide <- read.csv('olink_guide.csv')
# change subject # column nae
colnames(olink_guide)[1] <- "subjectnumber"
# select only 1y samples
olink_guide_1y <- subset(olink_guide, time == '1Y')
olink_data_1y <- merge(olink_data, olink_guide_1y, by = 'SampleID')
head(olink_data_1y)
subjectnumber_to_ID$subjectnumber <- subjectnumber_to_ID$Participant_ID
my_olink_data <- merge(olink_data_1y, subjectnumber_to_ID, by = 'subjectnumber')
my_olink_data <- subset(my_olink_data, select = -c(99:104))
my_olink_data
my_olink_data <- subset(my_olink_data, select = -c(1:2))
my_olink_data
row.names(my_olink_data) <- my_olink_data$CHILD_ID
my_olink_data <- subset(my_olink_data, select = -97)
my_olink_data <- subset(my_olink_data, select = -97)
my_olink_data

# now prepare olink data for maaslin2
head(my_3m_species)
my_olink_data_clean <- my_olink_data %>%
  select(-contains("Ctrl"))
head(my_olink_data_clean)
my_olink_data_clean$CHILD_ID <- row.names(my_olink_data_clean)
head(sg_my_metadata_3m)
head(my_olink_data_clean)
sample_id_child_id <- subset(sg_my_metadata_3m, select = c('SampleID', 'CHILD_ID'))
head(sample_id_child_id)
my_olink_data_clean <- merge(my_olink_data_clean, sample_id_child_id, by = 'CHILD_ID')
row.names(my_olink_data_clean) <- my_olink_data_clean$SampleID
my_olink_data_clean <- subset(my_olink_data_clean, select = -94)
my_olink_data_clean <- subset(my_olink_data_clean, select = -1)
olink_t <- t(my_olink_data_clean)
olink_t <- as.data.frame(olink_t)
class(my_3m_species$X7.048205.1.4)
class(olink_t$X7.048205.1.4)

#run maaslin2 on olink data and total IgA in milk
fit_data = Maaslin2(
  input_data = olink_t, 
  input_metadata = all_3m_metadata,
  output = "output_olink_iga",
  max_significance = 0.2,
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("group"),
  fixed_effects =  c("total_iga", "group") )

olink_total_iga_results_all <- read_tsv('output_olink_iga/all_results.tsv')
olink_total_iga_results <- subset(olink_total_iga_results_all, metadata == 'total_iga')
olink_total_iga_results <- subset(olink_total_iga_results, pval < 0.05)
#filtering for significant results
#plot Maaslin2 for significantly associated cytokines
olink_total_iga_plot <- plot_maaslin(olink_total_iga_results)
ggsave('olink_total_iga.svg', plot = olink_total_iga_plot)
# THIS IS FIGURE S1D


### RNA-Seq analysis in R
# RNAseq analysis May 2024
setwd("~/OneDrive - UBC/Documents/ECR characterization/RNAseq-fastqs")

library(tidyverse)
library(gridExtra)
library(ggplot2)
library(DESeq2)
library(rlang)
library(EnhancedVolcano)
library(textshaping)
library(gridExtra)
library(gtable)

#First, Specify the comparison you would like to look at. 
# (reference) vs (treatment)
# 1) CON vs ECR
# 2) CON vs HK
# 3) HK vs 

# read all read count tables (reads for each gene in each sample - 4 samples per group)
HK4tab <- read.table('HK4ReadsPerGene.out.tab')
write.csv(HK4tab, 'HK4tab.csv')  

HK3tab <- read.table('HK3ReadsPerGene.out.tab')
write.csv(HK3tab, 'HK3tab.csv')  

HK2tab <- read.table('HK2ReadsPerGene.out.tab')
write.csv(HK2tab, 'HK2tab.csv') 

HK1tab <- read.table('HK1ReadsPerGene.out.tab')
write.csv(HK1tab, 'HK1tab.csv') 

ECR1tab <- read.table('ECR1ReadsPerGene.out.tab')
write.csv(ECR1tab, 'ECR1tab.csv')  

ECR2tab <- read.table('ECR2ReadsPerGene.out.tab')
write.csv(ECR2tab, 'ECR2tab.csv')  

ECR3tab <- read.table('ECR3ReadsPerGene.out.tab')
write.csv(ECR3tab, 'ECR3tab.csv')  

ECR4tab <- read.table('ECR4ReadsPerGene.out.tab')
write.csv(ECR4tab, 'ECR4tab.csv')  

CON1tab <- read.table('CON1ReadsPerGene.out.tab')
write.csv(CON1tab, 'CON1tab.csv')  

CON2tab <- read.table('CON2ReadsPerGene.out.tab')
write.csv(CON2tab, 'CON2tab.csv')

CON3tab <- read.table('CON3ReadsPerGene.out.tab')
write.csv(CON3tab, 'CON3tab.csv')

CON4tab <- read.table('CON4ReadsPerGene.out.tab')
write.csv(CON4tab, 'CON4tab.csv')

# remove V3 and V4 (which are just repeats of V2-gene count)
HK1tab <- subset(HK1tab, select = -c(V3, V4))
colnames(HK1tab) <- c('gene_ID', 'HK1')
HK1tab <- HK1tab[-c(1:4), ]
HK2tab <- subset(HK2tab, select = -c(V3, V4))
colnames(HK2tab) <- c('gene_ID', 'HK2')
HK2tab <- HK2tab[-c(1:4), ]

HK3tab <- subset(HK3tab, select = -c(V3, V4))
colnames(HK3tab) <- c('gene_ID', 'HK3')
HK3tab <- HK3tab[-c(1:4), ]

HK4tab <- subset(HK4tab, select = -c(V3, V4))
colnames(HK4tab) <- c('gene_ID', 'HK4')
HK4tab <- HK4tab[-c(1:4), ]

ECR1tab <- subset(ECR1tab, select = -c(V3, V4))
colnames(ECR1tab) <- c('gene_ID', 'ECR1')
ECR1tab <- ECR1tab[-c(1:4), ]

ECR2tab <- subset(ECR2tab, select = -c(V3, V4))
colnames(ECR2tab) <- c('gene_ID', 'ECR2')
ECR2tab <- ECR2tab[-c(1:4), ]

ECR3tab <- subset(ECR3tab, select = -c(V3, V4))
colnames(ECR3tab) <- c('gene_ID', 'ECR3')
ECR3tab <- ECR3tab[-c(1:4), ]

ECR4tab <- subset(ECR4tab, select = -c(V3, V4))
colnames(ECR4tab) <- c('gene_ID', 'ECR4')
ECR4tab <- ECR4tab[-c(1:4), ]

CON1tab <- subset(CON1tab, select = -c(V3, V4))
colnames(CON1tab) <- c('gene_ID', 'CON1')
CON1tab <- CON1tab[-c(1:4), ]

CON2tab <- subset(CON2tab, select = -c(V3, V4))
colnames(CON2tab) <- c('gene_ID', 'CON2')
CON2tab <- CON2tab[-c(1:4), ]

CON3tab <- subset(CON3tab, select = -c(V3, V4))
colnames(CON3tab) <- c('gene_ID', 'CON3')
CON3tab <- CON3tab[-c(1:4), ]

CON4tab <- subset(CON4tab, select = -c(V3, V4))
colnames(CON4tab) <- c('gene_ID', 'CON4')
CON4tab <- CON4tab[-c(1:4), ]

gene_counts <- merge(ECR1tab, ECR2tab, by = "gene_ID", all = TRUE)
gene_counts <- merge(gene_counts, ECR3tab, by = "gene_ID")
gene_counts <- merge(gene_counts, ECR4tab, by = "gene_ID")
gene_counts <- merge(gene_counts, CON1tab, by = "gene_ID")
gene_counts <- merge(gene_counts, CON2tab, by = "gene_ID")
gene_counts <- merge(gene_counts, CON3tab, by = "gene_ID")
gene_counts <- merge(gene_counts, CON4tab, by = "gene_ID")
gene_counts <- merge(gene_counts, HK1tab, by = "gene_ID")
gene_counts <- merge(gene_counts, HK2tab, by = "gene_ID")
gene_counts <- merge(gene_counts, HK3tab, by = "gene_ID")
gene_counts <- merge(gene_counts, HK4tab, by = "gene_ID")
head(gene_counts)
# now I have a table for each sample & gene count for each gene

rownames(gene_counts) <- gene_counts$gene_ID
head(gene_counts)
gene_counts = gene_counts[,-1]
write.csv(gene_counts, 'gene_counts.csv')

# now I want to run DESeq2 for each comparison

# specify comparison
comparison = '1'

# define comparisons
if  (comparison == 1){
  #filtering gene_counts for desired samples (controls and ECR)
  GeneCounts =  gene_counts[,c(5:8,1:4)]
  
  #Setting metadata variables
  sample_names = colnames(GeneCounts)
  treatment_var = c("con","con","con","con","ECR","ECR", "ECR","ECR")
  
  #Setting variables for DESeq command
  relevel_factor = "con"
  design_factor = "treatment"
  comparison ="convsECR"
  
} else if (comparison == 2){
  #filtering the AllCounts database for desired samples===
  GeneCounts =  gene_counts[,c(5:12)]
  
  #Setting metadata variables
  sample_names = colnames(GeneCounts)
  treatment_var = c("con","con","con","con","HK","HK","HK","HK")
  
  #Setting variables for DESeq command
  relevel_factor = "con"
  design_factor = "treatment"
  comparison ="convsHK"
  
} else if (comparison == 3){ 
  #filtering the AllCounts database for desired samples
  GeneCounts =  gene_counts[,c(9:12,1:4)]
  
  #Setting metadata variables
  sample_names = colnames(GeneCounts)
  treatment_var = c("HK","HK","HK","HK","ECR","ECR","ECR","ECR")
  
  #Setting variables for DESeq command
  relevel_factor = "HK"
  design_factor = "treatment"
  comparison ="HKvsECR"
  
} 
# set metadata for DESeq2
metadata = data.frame(sample = sample_names,
                      treatment = treatment_var)
rownames(metadata) = metadata$sample


head(metadata)
dir.create(paste(comparison))
dds <- DESeqDataSetFromMatrix(countData = GeneCounts,
                              colData = metadata,
                              design = as.formula(paste("~",design_factor)))
#remove counts lower than 1
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]
#reference factor
dds[[design_factor]] <- relevel(dds[[design_factor]],relevel_factor)
dds_DESEQ = DESeq(dds)
dds

setwd(paste(comparison))

#Generate the expression file
dds_es <- estimateSizeFactors(dds)
expression_values = counts(dds_es, normalized=TRUE)
expression_values = cbind(gene_ID=rownames(expression_values),desc="na",expression_values)
colnames(expression_values)[1:2] = c("NAME","DESCRIPTION")

write.table(expression_values, file=paste(comparison,"_expression_values.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.csv(metadata,file = paste(comparison,"_metadata.csv"),
          row.names = FALSE, quote = FALSE)


#getting the raw results of differential expression
res <- results(dds_DESEQ)
res_example = as.data.frame(res)
res_PNA <- as.data.frame(res[,c(2:5)])
write.csv(res_PNA, file=paste(comparison,"_Raw_results.csv"))
mcols(res, use.names=TRUE)

#basic stats
summary(res)

#checking how many genes
length(res$padj)
res <- as.data.frame(res)
#how many sig genes
sum(res$padj<0.05,na.rm=TRUE)
res_SIG = as.data.frame(res) %>%
  filter(pvalue < 0.05)
write.csv(res_SIG, file=paste(comparison, "_pvalue0.05_results.csv"))
res_SIG_LFC <- subset(res_SIG, abs(log2FoldChange) > 1.5)
head(res_SIG_LFC)
write.csv(res_SIG_LFC, 'res_SIG_LFC_CONvsECR.csv')
# used this csv file to look for immune related genes
metadata

coldata = metadata


#This code creates a volcano plot highlighting the differentially expressed 
# genes above and below a certain log2foldchange 


res.df = as.data.frame(res)
res.df$IDS = rownames(res.df)
head(res.df)
res.df.interest <- subset(res.df, IDS %in% genes_of_interest)
genes_of_interest2 <- c('CXCL3','PIGR', 'CXCL8','CCL2',
                        'CXCL2','FCAMR','CXCL1','CCL20', 'C3','FABP1',
                        'REG4', 'IL1A')
genes_of_interest2
# FABP1 and REG4 not showing up - cuz don't meet cutoff
select = c("")
res.df
volplot <- EnhancedVolcano(res.df,
                           x = "log2FoldChange",
                           y = "pvalue",
                           lab = ifelse(res.df$IDS %in% genes_of_interest2, res.df$IDS, ""),
                           pCutoff = 0.05, 
                           FCcutoff = 1.5,
                           caption = NULL,  # Removes the caption
                           pointSize = 1,
                           boxedLabels = TRUE,
                           drawConnectors = TRUE,
                           widthConnectors = 0.5,
                           labSize = 4.0,
                           max.overlaps = Inf,
                           ylim = c(0,50),
                           title = expression("Genes enriched by" ~ italic("E. ramosum")),    # Italicizes E. ramosum in the title
                           subtitle = NULL,       # Removes the EnhancedVolcano subtitle
                           legendPosition = "none")  # Removes the legend


plot(volplot)
ggsave('volcano_ECR_CON.svg', plot = volplot)
# THIS IS FIGURE 4A

# KEGG term analysis
# after putting ECR hits into KEGG database

kegg_enriched_ECR_CON <- read_tsv('enrichment.KEGG.tsv')
head(kegg_enriched_ECR_CON)
kegg_enriched_ECR_CON <- as.data.frame(kegg_enriched_ECR_CON)
head(kegg_enriched_ECR_CON)

# Step 2: Data pre-processing (if needed)
# For example, filter out unnecessary columns
kegg_ECR_CON <- kegg_enriched_ECR_CON[, c("term_description", "observed_gene_count", "strength")]
colnames(kegg_ECR_CON) <- c("term", "gene_count", "strength")


ggplot(kegg_ECR_CON, aes(x = term, y = gene_count)) +
  geom_bar(stat = "identity") +
  labs(title = "STRING Analysis Results",
       x = "Gene ID",
       y = "Interaction Score") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

kegg_ECR_CON

head(kegg_enriched_ECR_CON)
colnames(kegg_enriched_ECR_CON) <- c('term_id', 'term_description', 'observed_gene_count', 'background_gene_count', 
                                     'strength', 'false_discovery_rate', 'matching_proteins', 'matching_protein_labels')
class(kegg_enriched_ECR_CON$matching_proteins)


head(kegg_enriched_ECR_CON)
kegg_enriched_ECR_CON$pathway_weight = kegg_enriched_ECR_CON$background_gene_count/kegg_enriched_ECR_CON$observed_gene_count

kegg_enriched_ECR_CON$term_description = factor(kegg_enriched_ECR_CON$term_description, levels = kegg_enriched_ECR_CON[order(kegg_enriched_ECR_CON$strength),"term_description"])

pathway_plot_kegg = ggplot(data=kegg_enriched_ECR_CON,aes(x = strength,y=term_description, color = false_discovery_rate ))+
  geom_point(aes(size = observed_gene_count))+
  theme_classic()+
  scale_color_gradient(low = "cornflowerblue",high = "brown1") +
  labs(y = "KEGG pathway", size = "observed gene count", color = 'FDR')
pathway_plot_kegg
ggsave(filename ="pathway_plot_kegg.tiff" ,plot = pathway_plot_kegg, dpi = 600, width = 8, height = 6)
# THIS IS FIGURE 4B


rownames(all_3m_metadata) <- all_3m_metadata$SampleID
head(my_3m_species)
mean(colSums(my_3m_species))
# subsetting species data to abundance cut off of 0.001 (0.1% RA)
Blongum_infants <- my_3m_species[, 
                my_3m_species["Bifidobacterium_longum", ] > 0 ]
                
fit_data = Maaslin2(
  input_data = Blongum_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_blongum_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("binf_iga", "Exact_age_months"))

ssalivarius_infants <- my_3m_species[, 
                                 my_3m_species["Streptococcus_parasanguinis", ] > 0 ]

fit_data = Maaslin2(
  input_data = ssalivarius_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_ssalivarius_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("ssalivarius_iga", "Exact_age_months"))

sparasan_infants <- my_3m_species[, 
                                     my_3m_species["Streptococcus_parasanguinis", ] > 0 ]

fit_data = Maaslin2(
  input_data = sparasan_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_sparasan_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("sparasan_iga", "Exact_age_months"))

rumino_infants <- my_3m_species[, 
                                my_3m_species["Ruminococcus_gnavus", ] > 0 ] 
fit_data = Maaslin2(
  input_data = rumino_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_rumino_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("rumino_iga", "Exact_age_months"))


para_infants <- my_3m_species[, 
                                my_3m_species["Parabacteroides_distasonis", ] > 0 ]
fit_data = Maaslin2(
  input_data = para_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_para_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("para_iga", "Exact_age_months"))

head(all_3m_metadata)
kvariicola_infants <- my_3m_species[,
                  my_3m_species["Klebsiella_variicola", ] > 0 ]

fit_data = Maaslin2(
  input_data = kvariicola_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_kvariicola_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("kvari_iga", "Exact_age_months"))
                  

ecoli_infants <- my_3m_species[,
                                    my_3m_species["Escherichia_coli", ] > 0 ]

fit_data = Maaslin2(
  input_data = ecoli_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_ecoli_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("ecoli_iga", "Exact_age_months"))

                  
eramosum_infants <- my_3m_species[, 
                               my_3m_species["Erysipelatoclostridium_ramosum", ] > 0 ]
fit_data = Maaslin2(
  input_data = eramosum_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_eramosum_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("eramosum_iga", "Exact_age_months"))


                  
entero_infants <- my_3m_species[,
                                  my_3m_species["Enterococcus_faecalis", ] > 0 ]
                  
fit_data = Maaslin2(
  input_data = entero_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_entero_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("entero_iga", "Exact_age_months"))

kpn_infants <- my_3m_species[,
                                my_3m_species["Klebsiella_pneumoniae", ] > 0 ]

fit_data = Maaslin2(
  input_data = kpn_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_kpn_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("kpn_iga", "Exact_age_months"))



cperfing_infants <- my_3m_species[,
                                my_3m_species["Clostridium_perfringens", ] > 0 ]


fit_data = Maaslin2(
  input_data = cperfing_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_cperfing_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("cperfing_iga", "Exact_age_months"))


               
bbrev_infants <- my_3m_species[,
                                  my_3m_species["Bifidobacterium_breve", ] > 0 ]

fit_data = Maaslin2(
  input_data = bbrev_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_bbrev_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("bbrev_iga", "Exact_age_months"))

                  
bbifidum_infants <- my_3m_species[,
                               my_3m_species["Bifidobacterium_bifidum", ] > 0 ]
fit_data = Maaslin2(
  input_data = bbifidum_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_bbifidum_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("bbifidum_iga", "Exact_age_months"))


                  
buniform_infants <- my_3m_species[,
                                  my_3m_species["Bacteroides_uniformis", ] > 0 ]
fit_data = Maaslin2(
  input_data = buniform_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_buniform_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("buniform_iga", "Exact_age_months"))


bfrag_infants <- my_3m_species[,
                                  my_3m_species["Bacteroides_fragilis", ] > 0 ]
fit_data = Maaslin2(
  input_data = bfrag_infants, 
  input_metadata = all_3m_metadata,
  output = "output_strat_bfrag_iga",
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  random_effects = ("site"),
  fixed_effects = c("bfrag_iga", "Exact_age_months"))


binf_iga_results_strat <- read_tsv('output_strat_blongum_iga/all_results.tsv')
binf_iga_results_strat <- subset(binf_iga_results_strat, metadata == 'binf_iga')
bbrev_iga_results_strat <- read_tsv('output_strat_bbrev_iga/all_results.tsv')
bbrev_iga_results_strat<- subset(bbrev_iga_results_strat, metadata == 'bbrev_iga')
bbif_iga_results_strat <- read_tsv('output_strat_bbifidum_iga/all_results.tsv')
bbif_iga_results_strat <- subset(bbif_iga_results_strat, metadata == 'bbifidum_iga')
ecoli_iga_results_strat <- read_tsv('output_strat_ecoli_iga/all_results.tsv')
ecoli_iga_results_strat <- subset(ecoli_iga_results_strat, metadata == 'ecoli_iga')
kvari_iga_results_strat <- read_tsv('output_strat_kvariicola_iga/all_results.tsv')
kvari_iga_results_strat <- subset(kvari_iga_results_strat, metadata == 'kvari_iga')
kpn_iga_results_strat <- read_tsv('output_strat_kpn_iga/all_results.tsv')
kpn_iga_results_strat <- subset(kpn_iga_results_strat, metadata == 'kpn_iga')
bfrag_iga_results_strat<- read_tsv('output_strat_bfrag_iga/all_results.tsv')
bfrag_iga_results_strat<- subset(bfrag_iga_results_strat, metadata == 'bfrag_iga')
buniform_iga_results_strat <- read_tsv('output_strat_buniform_iga/all_results.tsv')
buniform_iga_results_strat <- subset(buniform_iga_results_strat, metadata == 'buniform_iga')
rumino_iga_results_strat<- read_tsv('output_strat_rumino_iga/all_results.tsv')
rumino_iga_results_strat<- subset(rumino_iga_results_strat, metadata == 'rumino_iga')
entero_iga_results_strat<- read_tsv('output_strat_entero_iga/all_results.tsv')
entero_iga_results_strat<- subset(entero_iga_results_strat, metadata == 'entero_iga')
eramosum_iga_results_strat<- read_tsv('output_strat_eramosum_iga/all_results.tsv')
eramosum_iga_results_strat<- subset(eramosum_iga_results_strat, metadata == 'eramosum_iga')
cperfing_iga_results_strat<- read_tsv('output_strat_cperfing_iga/all_results.tsv')
cperfing_iga_results_strat<- subset(cperfing_iga_results_strat, metadata == 'cperfing_iga')
ssalivarius_iga_results_strat<- read_tsv('output_strat_ssalivarius_iga/all_results.tsv')
ssalivarius_iga_results_strat<- subset(ssalivarius_iga_results_strat, metadata == 'ssalivarius_iga')
sparasan_iga_results_strat<- read_tsv('output_strat_sparasan_iga/all_results.tsv')
sparasan_iga_results_strat<- subset(sparasan_iga_results_strat, metadata == 'sparasan_iga')
para_iga_results_strat<- read_tsv('output_strat_para_iga/all_results.tsv')
para_iga_results_strat<- subset(para_iga_results_strat, metadata == 'para_iga')

# individuate each data frame so that all can be merged
strat_species_binf <- subset(binf_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_binf) <- c('feature', 'binf_coef', 'binf_pval', 'binf_stderr')
strat_species_bbif <- subset(bbif_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_bbif) <- c('feature', 'bbif_coef', 'bbif_pval', 'bbif_stderr')
strat_species_bbrev <- subset(bbrev_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_bbrev) <- c('feature', 'bbrev_coef', 'bbrev_pval', 'bbrev_stderr')
strat_species_ecoli <- subset(ecoli_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_ecoli) <- c('feature', 'ecoli_coef', 'ecoli_pval', 'ecoli_stderr')
strat_species_kpn <- subset(kpn_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_kpn) <- c('feature', 'kpn_coef', 'kpn_pval', 'kpn_stderr')
strat_species_kvari <- subset(kvari_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_kvari) <- c('feature', 'kvari_coef', 'kvari_pval', 'kvari_stderr')
strat_species_bfrag <- subset(bfrag_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_bfrag) <- c('feature', 'bfrag_coef', 'bfrag_pval', 'bfrag_stderr')
strat_species_buniform <- subset(buniform_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_buniform) <- c('feature', 'buniform_coef', 'buniform_pval', 'buniform_stderr')
strat_species_rumino <- subset(rumino_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_rumino) <- c('feature', 'rumino_coef', 'rumino_pval', 'rumino_stderr')
strat_species_entero <- subset(entero_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_entero) <- c('feature', 'entero_coef', 'entero_pval', 'entero_stderr')
strat_species_eramosum <- subset(eramosum_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_eramosum) <- c('feature', 'eramosum_coef', 'eramosum_pval', 'eramosum_stderr')
strat_species_cperfing <- subset(cperfing_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_cperfing) <- c('feature', 'cperfing_coef', 'cperfing_pval', 'cperfing_stderr')
strat_species_ssalivarius <- subset(ssalivarius_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_ssalivarius) <- c('feature', 'ssalivarius_coef', 'ssalivarius_pval', 'ssalivarius_stderr')
strat_species_sparasan <- subset(sparasan_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_sparasan) <- c('feature', 'sparasan_coef', 'sparasan_pval', 'sparasan_stderr')
strat_species_para <- subset(para_iga_results_strat, select = c('feature', 'coef', 'pval', 'stderr'))
colnames(strat_species_para) <- c('feature', 'para_coef', 'para_pval', 'para_stderr')



strat_species_table <- merge(strat_species_bfrag, strat_species_bbrev, by = 'feature')
head(strat_species_table)
strat_species_table <- merge(strat_species_table, strat_species_binf, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_bbif, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_buniform, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_ecoli, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_kpn, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_kvari, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_sparasan, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_ssalivarius, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_rumino, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_entero, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_eramosum, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_cperfing, by = 'feature')
strat_species_table <- merge(strat_species_table, strat_species_para, by = 'feature')

strat_species_table_coef <- subset(strat_species_table, select = c(feature, binf_coef, bbrev_coef, bbif_coef, 
                                                       bfrag_coef, buniform_coef, ecoli_coef, kpn_coef, 
                                                       kvari_coef, sparasan_coef, ssalivarius_coef, rumino_coef, 
                                                       entero_coef, eramosum_coef, cperfing_coef, para_coef))
colnames(strat_species_table_coef) <- c('feature',  'B. longum IgA', 'B. breve IgA', 'B. bifidum IgA', 
                                  'B. fragilis IgA', 'B. uniformis IgA', 'E. coli IgA', 'K. pneumoniaeae IgA', 
                                  'K. variicola IgA', 'S. parasanguinis IgA', 'S. salivarius IgA', 
                                  'R. gnavus IgA', 'E. faecalis IgA', 'E. ramosum IgA', 'C. perfingens IgA',
                                  'P. disasonis IgA')
head(strat_species_table_coef)
strat_data_melt <- melt(strat_species_table_coef)
head(strat_data_melt)


# now want to get all species-specific IgA associations and total IgA associations for a bubble plot (Figure 1C)
colnames(strat_species_bbrev) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igabbrev <- subset(strat_species_bbrev, feature == 'Bifidobacterium_breve')
strat_igabbrev <- as.data.frame(strat_igabbrev)
colnames(strat_species_bbif) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igabbifidum <- subset(strat_species_bbif, feature == 'Bifidobacterium_bifidum')
strat_igabbifidum <- as.data.frame(strat_igabbifidum)
colnames(strat_species_binf) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igabinf <- subset(strat_species_binf, feature == 'Bifidobacterium_longum')
strat_igabinf <- as.data.frame(strat_igabinf)
colnames(strat_species_ecoli) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igaecoli <- subset(strat_species_ecoli, feature == 'Escherichia_coli')
strat_igaecoli <- as.data.frame(strat_igaecoli)
colnames(strat_species_kpn) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igakpn <- subset(strat_species_kpn, feature == 'Klebsiella_pneumoniae')
strat_igakpn <- as.data.frame(strat_igakpn)
colnames(strat_species_kvari) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igakvari <- subset(strat_species_kvari, feature == 'Klebsiella_variicola')
strat_igakvari <- as.data.frame(strat_igakvari)
colnames(strat_species_bfrag) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igabfrag <- subset(strat_species_bfrag, feature == 'Bacteroides_fragilis')
strat_igabfrag <- as.data.frame(strat_igabfrag)
colnames(strat_species_buniform) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igabuniform <- subset(strat_species_buniform, feature == 'Bacteroides_uniformis')
strat_igabuniform <- as.data.frame(strat_igabuniform)
colnames(strat_species_rumino) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igarumino <- subset(strat_species_rumino, feature == 'Ruminococcus_gnavus')
strat_igarumino <- as.data.frame(strat_igarumino)
colnames(strat_species_entero) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igaentero <- subset(strat_species_entero, feature == 'Enterococcus_faecalis')
strat_igaentero <- as.data.frame(strat_igaentero)
colnames(strat_species_eramosum) <- c('feature', 'spec_coef', 'spec_pval', 'spec_stderr')
strat_igaeramosum <- subset(strat_species_eramosum, feature == 'Erysipelatoclostridium_ramosum')
strat_igaeramosum <- as.data.frame(strat_igaeramosum)
head(species_cperfing)
colnames(species_cperfing) <- c('feature', 'spec_coef', 'spec_stderr')
strat_stacked_df_con <- rbind(strat_igabbifidum, strat_igabbrev,
                              strat_igabfrag, strat_igabinf, strat_igabuniform,
                              strat_igaecoli, strat_igaentero, strat_igaeramosum,
                              strat_igakpn, strat_igakvari, strat_igarumino)
strat_stacked_df_con <- strat_stacked_df_con %>% mutate(spec_sig = ifelse(abs(spec_pval) < 0.1, 'y', 'n'))
# none were significant doing this - probably cuz IgA has an effect on presence absence as 
# well as relative abundance
stacked_df_con

              
# no infant had all of the species


### R code for Figure 2B and Figure 2C using data on the entire CHILD cohort:

load("Microbiota_Metadata_Kate.RData") # includes metadata, taxa prevalence data
ramosum<-taxa.table.prev[grep("ramosum", rownames(taxa.table.prev)), ] %>% t()
colnames(ramosum)<-gsub(" .*", "",colnames(ramosum))
colnames(ramosum)
metadata.rbind<-merge(metadata.rbind, ramosum, by="row.names")
rownames(metadata.rbind)<-metadata.rbind$`#SampleID`

metadata.rbind$allergy.case<-as.factor(metadata.rbind$allergy.case)

#Creating a psuedocount that is one half the minimumum of E.ramosum. 
#This will only be used to create log10 graphs for easier visuals
pseudo<-(min(metadata.rbind$Erysipelatoclostridium_ramosum[!metadata.rbind$Erysipelatoclostridium_ramosum==0])/2)


## Figue 2B: comparison of E. ramosum at 3mo and 1y by allergy status

# create an allergy label
metadata.rbind$allergy.label <- ifelse(metadata.rbind$allergy.case == 1, "allergic", "control")

# Create an interaction between age and allergy status using the new labels
metadata.rbind$Visit_allergy <- factor(
  interaction(metadata.rbind$Visit, metadata.rbind$allergy.label, sep = "."),
  levels = c("3 month.control", "1 year.control", "3 month.allergic", "1 year.allergic")
)

my_comparisons <- list(
  c("3 month.control", "1 year.control"),
  c("3 month.allergic", "1 year.allergic"),
  c("3 month.control", "3 month.allergic"), 
  c("1 year.control", "1 year.allergic")
)
# Then use it in your ggplot
over_time_plot <- ggplot(
  subset(metadata.rbind, !is.na(allergy.case) & bf_3m_status %in% c("exclusive", "partial")),
  aes(x = Visit_allergy, y = Erysipelatoclostridium_ramosum + pseudo, fill = allergy.case)
) +
  geom_boxplot(color = "black") +
  geom_point(position = position_jitter(width = 0.2), color = "black", alpha = 0.6) +
  scale_y_continuous(
    trans = "log10",
    breaks = scales::log_breaks(base = 10),
    labels = scales::label_log(base = 10)
  ) +
  annotation_logticks(sides = "l") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme_classic() +
  theme(
    legend.position = "none",                    # remove legend
    axis.title.x = element_text(size = 16),      # increase x-axis label size
    axis.title.y = element_text(size = 16),      # increase y-axis label size
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  labs(
    x = "Visit, case/control",
    y = expression(italic("Erysipelatoclostridium ramosum")~" (log-transformed)"))
over_time_plot


# then use in violin plot:

over_time_plot_violin <- ggplot(
  subset(metadata.rbind, !is.na(allergy.case) & bf_3m_status %in% c("exclusive", "partial")),
  aes(x = Visit_allergy, y = Erysipelatoclostridium_ramosum + pseudo, fill = allergy.case)
) +
  geom_violin(color = "black") +
  geom_point(position = position_jitter(width = 0.2), color = "black", alpha = 0.6) +
  scale_y_continuous(
    trans = "log10",
    breaks = scales::log_breaks(base = 10),
    labels = scales::label_log(base = 10)
  ) +
  annotation_logticks(sides = "l") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme_classic() +
  theme(
    legend.position = "none",                    # remove legend
    axis.title.x = element_text(size = 16),      # increase x-axis label size
    axis.title.y = element_text(size = 16),      # increase y-axis label size
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  labs(
    x = "Visit, case/control",
    y = expression(italic("Erysipelatoclostridium ramosum")~" (log-transformed)"))
over_time_plot_violin

ggsave('eramosum-over-time-violin.svg', plot = over_time_plot_violin) #Figure 2B


# Figure 2C: E. ramosum by breastfeeding status at 3 months
my_comparisons <- list( c(1, 2), c(1,3))
breastfeeding_eramosum <- ggplot(
  subset(metadata.rbind, !is.na(bf_3m_status) & Visit == "3 month"),
  aes(x = bf_3m_status, y = Erysipelatoclostridium_ramosum + pseudo, fill = bf_3m_status)
) +
  geom_boxplot(color = "black") +
  geom_point(position = position_jitter(width = 0.2), color = "black", alpha = 0.6) +
  scale_y_continuous(
    trans = "log10",
    breaks = scales::log_breaks(base = 10),
    labels = scales::label_log(base = 10)
  ) +
  annotation_logticks(sides = "l") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  labs(
    x = "breastfeeding status at 3 months",
    y = expression(italic("Erysipelatoclostridium ramosum")~"(log transformed)")
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "none"
  )
breastfeeding_eramosum

prevalence_bf <- metadata.rbind %>%
  filter(
    !is.na(bf_3m_status),
    Visit == "3 month"
  ) %>%
  group_by(bf_3m_status) %>%
  summarise(
    n_samples = n(),
    n_nonzero = sum(Erysipelatoclostridium_ramosum > 0, na.rm = TRUE),
    prevalence = n_nonzero / n_samples
  ) %>%
  ungroup()

ggsave('breastfeeding-eramosum.svg', plot = breastfeeding_eramosum) # Figure 2C