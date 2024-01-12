library(tidyverse)

setwd("/home/filip/Desktop/LOH/dip_and_trip_DP/R_DP_data")

transformed_ad_values <- read.csv("~/Desktop/LOH/dip_and_trip_DP/R_DP_data/transformed_ad_values.csv")
View(transformed_ad_values)

# instead samples names in many columns, transform them into only 2
ad_long <- transformed_ad_values %>% 
  pivot_longer(cols = everything(), 
               names_to = "names", 
               values_to = "values")


names(transformed_ad_values)

ggplot(transformed_ad_values, aes(x=csc102.vcf_0)) + 
  geom_histogram(binwidth = 1, # adjusting binwidth
                 fill="blue", 
                 color="black") + 
  theme_minimal() +
  labs(title="Histogram of Values", x="Value", y="Count")

# Individual density plot
ggplot(dp_values, aes(x=csc102.vcf)) + 
  geom_density(fill="blue", alpha=0.5) + 
  scale_x_continuous(limits = c(10, 1000)) +
  theme_minimal() +
  labs(title="Density Plot of Values", x="Value", y="Density")


# Multiple overlay density plot
# csc102
ad_long %>%
  filter(names %in% c('csc102.vcf_0', 'csc102.vcf_1')) %>%
  ggplot(aes(x = values, fill = names)) +
  geom_density(alpha = 0.5) +
  labs(x = 'AD', y = 'Density') +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 40))

# Multiple overlay density plot
# csc103
ad_long %>%
  filter(names %in% c('csc103.vcf_0', 'csc103.vcf_1')) %>%
  ggplot(aes(x = values, fill = names)) +
  geom_density(alpha = 0.5) +
  labs(x = 'AD', y = 'Density') +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 40))

# Multiple overlay density plot
# csc134
ad_long %>%
  filter(names %in% c('csc134.vcf_0', 'csc134.vcf_1')) %>%
  ggplot(aes(x = values, fill = names)) +
  geom_density(alpha = 0.5) +
  labs(x = 'AD', y = 'Density') +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 120))

# Multiple overlay density plot
# csc138
ad_long %>%
  filter(names %in% c('csc138.vcf_0', 'csc138.vcf_1')) %>%
  ggplot(aes(x = values, fill = names)) +
  geom_density(alpha = 0.5) +
  labs(x = 'AD', y = 'Density') +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 120))

# DIPLOIDS
# Multiple overlay density plot
# csc163
ad_long %>%
  filter(names %in% c('csc163.vcf_0', 'csc163.vcf_1')) %>%
  ggplot(aes(x = values, fill = names)) +
  geom_density(alpha = 0.5) +
  labs(x = 'AD', y = 'Density') +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 75))

# Multiple overlay density plot
# csc164
ad_long %>%
  filter(names %in% c('csc164.vcf_0', 'csc164.vcf_1')) %>%
  ggplot(aes(x = values, fill = names)) +
  geom_density(alpha = 0.5) +
  labs(x = 'AD', y = 'Density') +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 100))

# Multiple overlay density plot
# csc174
ad_long %>%
  filter(names %in% c('csc174.vcf_0', 'csc174.vcf_1')) %>%
  ggplot(aes(x = values, fill = names)) +
  geom_density(alpha = 0.5) +
  labs(x = 'AD', y = 'Density') +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 50))




ad_long %>%
  filter(names %in% c('csc164.vcf_0', 'csc164.vcf_1')) %>%
  group_by(names) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = names, values_from = values) %>%
  select(-row) %>%
  ggplot(aes(x = csc164.vcf_0, y = csc164.vcf_1)) +
  geom_point() + 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')


# TRANSFORM -> SUM AND THEN MEAN AND PLOT
csc164 <- ad_long %>%
  filter(names %in% c('csc164.vcf_0', 'csc164.vcf_1')) %>%
  group_by(names) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = names, values_from = values) %>%
  select(-row) 

# suma <- sum(csc164$csc164.vcf_0) + sum(csc164$csc164.vcf_1)
# csc164$csc164.vcf_0 <- csc164$csc164.vcf_0 / suma
# csc164$csc164.vcf_1 <- csc164$csc164.vcf_1 / suma




csc164$sum <- csc164$csc164.vcf_0 + csc164$csc164.vcf_1
csc164$csc164.vcf_0 <- csc164$csc164.vcf_0 / csc164$sum
csc164$csc164.vcf_1 <- csc164$csc164.vcf_1 / csc164$sum

csc164 %>%
  ggplot(aes(x = csc164.vcf_1, y = csc164.vcf_0))+
  geom_point()

csc164 %>%
  pivot_longer(cols = everything(), 
               names_to = "names", 
               values_to = "values") %>%
  filter(names %in% c('csc164.vcf_0', 'csc164.vcf_1')) %>%
  ggplot(aes(x = values, fill = names)) +
  geom_density(alpha = 0.5) +
  labs(title = "AD % values density plot",x = 'AD', y = 'Density') +
  theme_minimal() 
  


## ++ genotype k AD