library(ggplot2)

names(dp_values)

# Read the CSV file
df <- read.csv("/home/filip/Desktop/LOH/dip_and_trip_DP/R_DP_data/dp_values.csv")
dp_values <- read.csv("/home/filip/Desktop/LOH/dip_and_trip_DP/R_DP_data/dp_values.csv")

# Transform the data
long_df <- df %>% 
  pivot_longer(cols = everything(), 
               names_to = "names", 
               values_to = "values")


# Write the transformed data to a new CSV file
write.csv(long_df, "transformed_dp_values.csv", row.names = FALSE)


# Individual histogram
ggplot(dp_values, aes(x=csc102.vcf)) + 
  geom_histogram(binwidth = 1, # adjusting binwidth
                 fill="blue", 
                 color="black") + 
  scale_x_continuous(limits = c(10,60)) +
  theme_minimal() +
  labs(title="Histogram of Values", x="Value", y="Count")

# Individual density plot
ggplot(dp_values, aes(x=csc102.vcf)) + 
  geom_density(fill="blue", alpha=0.5) + 
  scale_x_continuous(limits = c(10, 1000)) +
  theme_minimal() +
  labs(title="Density Plot of Values", x="Value", y="Density")
  
  
# Multiple overlay density plot
  ggplot(long_df,
         aes(x = values,
             fill = names)) +
    geom_density(alpha = 0.9) +
    scale_x_continuous(limits = c(10,120)) +
    theme_minimal() +
    labs(x = "DP", y = "Density")

  

# filtr pozic z tabulky, které se navázali na vcfko hybridů 42 009 pozic
# vyfiltrovat název, pozici, referenci, alternativu
# když na např 0/1 a reference T alternativa G přelož a připočítej do tabulky
  # na danou pozici
  
  
  
  
  