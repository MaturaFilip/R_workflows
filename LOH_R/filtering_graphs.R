library(tidyverse)

extracted_ch_pos_refalt_sample_dp <- read.csv("~/Desktop/LOH/dip_and_trip_DP/R_DP_data/genotypy_pozice_tabulka/extracted_ch_pos_refalt_sample_dp.csv")
View(extracted_ch_pos_refalt_sample_dp)

extracted_ch_pos_refalt_sample_dp <- as_tibble(extracted_ch_pos_refalt_sample_dp)
typeof(extracted_ch_pos_refalt_sample_dp)


extracted_ch_pos_refalt_sample_dp<- extracted_ch_pos_refalt_sample_dp %>%
  filter(DP >= 10)


extracted_ch_pos_refalt_sample_dp <- extracted_ch_pos_refalt_sample_dp %>%
  select(-DP)
  


# Create a "new rows"
# split the genotype into separate rows (0/1/1) -> 0, 1, 1 in individual rows
df <- extracted_ch_pos_refalt_sample_dp %>%
  separate_longer_delim(Genotype, delim = '/')
# again but now with |
df <- df %>%
  separate_longer_delim(Genotype, delim = '|')




# second df with results for locus and position
# results <- extracted_ch_pos_refalt_sample_dp %>%
#   select(Chromosome, Position)

get_nucleotide <- function(genotype, reference, alternative) {
  if (genotype == 0) {
    return(reference)
  } else if (genotype == 1) {
    return(alternative)
  } else if (genotype == ".") {
    return("NONE")
  } else {
    return(NA) # for unexpected values
  }
}

# apply function to each row
df$Nucleotide <- mapply(get_nucleotide, df$Genotype, df$Reference, df$Alternative)


# if the alternative allele if longer than 1 delete the row
df <- df %>% 
  filter(nchar(Reference) == 1 & nchar(Alternative) == 1 & 
           Nucleotide %in% c("A", "C", "G", "T") | Reference == ".")

df$set3_A <- 0
df$set3_C <- 0
df$set3_G <- 0
df$set3_T <- 0


# Modify column A,C,G,T based on Nucleotide column
# A -> if Nucleotide == A then return 1, if not, return 0
df <- df %>%
  mutate(set3_A = set3_A + case_when(Nucleotide == "A" ~ 1, TRUE ~ 0),
         set3_C = set3_C + case_when(Nucleotide == "C" ~ 1, TRUE ~ 0),
         set3_G = set3_G + case_when(Nucleotide == "G" ~ 1, TRUE ~ 0),
         set3_T = set3_T + case_when(Nucleotide == "T" ~ 1, TRUE ~ 0))



# split the data into groups based on one or more variable
# summarise -> reduce the rows into one and apply function (here sum())
df_all<- df %>%
  group_by(Chromosome, Position) %>%
  summarise(set3_A = sum(set3_A), 
            set3_C = sum(set3_C),
            set3_G = sum(set3_G),
            set3_T = sum(set3_T))


write.csv(df, "/home/filip/Desktop/LOH/dip_and_trip_DP/R_DP_data/genotypy_pozice_tabulka/genotypes_sum_all.csv", row.names=FALSE)



##############################################################################
##############################################################################
##############################################################################
##################SAME, BUT FOR DIPLOIDS######################################
##############################################################################
##############################################################################

df_diploids <- extracted_ch_pos_refalt_sample_dp %>%
  filter(Sample == "csc163" | Sample == "csc164" | Sample == "csc174")


df_diploids <- df_diploids %>%
  separate_longer_delim(Genotype, delim = '/')
# again but now with |
df_diploids<- df_diploids %>%
  separate_longer_delim(Genotype, delim = '|')


df_diploids$Nucleotide <- mapply(get_nucleotide, df_diploids$Genotype,
                                 df_diploids$Reference,
                                 df_diploids$Alternative)



df_diploids <- df_diploids %>% 
  filter(nchar(Reference) == 1 & nchar(Alternative) == 1 & 
           Nucleotide %in% c("A", "C", "G", "T") | Reference == ".")



df_diploids$set3_A <- 0
df_diploids$set3_C <- 0
df_diploids$set3_G <- 0
df_diploids$set3_T <- 0


# Modify column A,C,G,T based on Nucleotide column
# A -> if Nucleotide == A then return 1, if not, return 0
df_diploids <- df_diploids %>%
  mutate(set3_A = set3_A + case_when(Nucleotide == "A" ~ 1, TRUE ~ 0),
         set3_C = set3_C + case_when(Nucleotide == "C" ~ 1, TRUE ~ 0),
         set3_G = set3_G + case_when(Nucleotide == "G" ~ 1, TRUE ~ 0),
         set3_T = set3_T + case_when(Nucleotide == "T" ~ 1, TRUE ~ 0))


# split the data into groups based on one or more variable
# summarise -> reduce the rows into one and apply function (here sum())
df_diploids <- df_diploids %>%
  group_by(Chromosome, Position) %>%
  summarise(set3_A = sum(set3_A), 
            set3_C = sum(set3_C),
            set3_G = sum(set3_G),
            set3_T = sum(set3_T))


write.csv(df_diploids, "/home/filip/Desktop/LOH/dip_and_trip_DP/R_DP_data/genotypy_pozice_tabulka/genotypes_sum_diploids.csv", row.names=FALSE)


##############################################################################
##############################################################################
##############################################################################
##################SAME, BUT FOR TRIPLOIDS#####################################
##############################################################################
##############################################################################


df_triploids <- extracted_ch_pos_refalt_sample_dp %>%
  filter(Sample == "csc102" | Sample == "csc103" | Sample == "csc134" | Sample == "csc138")

df_triploids <- df_triploids %>%
  separate_longer_delim(Genotype, delim = '/')
# again but now with |
df_triploids<- df_triploids %>%
  separate_longer_delim(Genotype, delim = '|')


df_triploids$Nucleotide <- mapply(get_nucleotide, df_triploids$Genotype,
                                  df_triploids$Reference,
                                  df_triploids$Alternative)



df_triploids <- df_triploids %>% 
  filter(nchar(Reference) == 1 & nchar(Alternative) == 1 & 
           Nucleotide %in% c("A", "C", "G", "T") | Reference == ".")



df_triploids$set3_A <- 0
df_triploids$set3_C <- 0
df_triploids$set3_G <- 0
df_triploids$set3_T <- 0


# Modify column A,C,G,T based on Nucleotide column
# A -> if Nucleotide == A then return 1, if not, return 0
df_triploids <- df_triploids %>%
  mutate(set3_A = set3_A + case_when(Nucleotide == "A" ~ 1, TRUE ~ 0),
         set3_C = set3_C + case_when(Nucleotide == "C" ~ 1, TRUE ~ 0),
         set3_G = set3_G + case_when(Nucleotide == "G" ~ 1, TRUE ~ 0),
         set3_T = set3_T + case_when(Nucleotide == "T" ~ 1, TRUE ~ 0))


# split the data into groups based on one or more variable
# summarise -> reduce the rows into one and apply function (here sum())
df_triploids <- df_triploids %>%
  group_by(Chromosome, Position) %>%
  summarise(set3_A = sum(set3_A), 
            set3_C = sum(set3_C),
            set3_G = sum(set3_G),
            set3_T = sum(set3_T))


write.csv(df_triploids, "/home/filip/Desktop/LOH/dip_and_trip_DP/R_DP_data/genotypy_pozice_tabulka/genotypes_sum_triploids.csv", row.names=FALSE)










