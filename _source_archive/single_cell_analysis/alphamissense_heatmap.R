#clinvar and alphamissense heatmap 

salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)





data_split <- data %>%
  mutate(
    Reference_amino_acid = substr(V2, 1, 1),
    mutated_amino_acid = substr(V2, nchar(V2), nchar(V2)),
    V2 = substr(V2, 2, nchar(V2) - 1)  # Keep the middle part of V2
  )
# View the result
head(data_split)




      # Combine Reference_amino_acid and V2 into a new column
data_split <- data_split %>%
  mutate(Reference_V2 = paste(Reference_amino_acid, V2, sep = " "))

# View the first few rows of the updated dataframe
head(data_split)



data_split$Reference_V2 <- factor(data_split$Reference_V2, levels = unique(data_split$Reference_V2))




library(ggplot2)

# Create a heatmap with ggplot
ggplot(data_split, aes(x = Reference_V2, y = mutated_amino_acid, fill = V3)) +
  geom_tile(colour="black") +
  scale_fill_gradientn(colours=c('blue', 'white', 'red')) +  # No white, direct transition from blue to red
  labs(title = "Heatmap of Reference vs Mutated Amino Acids HOPX Alphamissense", 
       x = "Reference Amino Acid", 
       y = "Mutated Amino Acid", 
       fill = "Score (V3)") + 
  #scale_x_discrete(labels = x) +  # Custom labels for x-axis
  theme_minimal() +  # Use a minimal theme for a cleaner look
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background=element_rect(fill='black'), panel.grid=element_blank())  # Rotate x-axis labels for better readability




salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=10GB --job-name=TEST --time=00:15:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)


library(ggplot2)
library(dplyr)

# Create a heatmap with ggplot
ggplot(merged, aes(x = Reference_V2, y = mutated_amino_acid, fill = score)) +  # Ensure V3 is the pathogenic score
  geom_tile(colour = "black") +  # Create tiles for the heatmap
  scale_fill_gradientn(colours = c('blue', 'white', 'red')) +  # Color gradient
  labs(title = "Heatmap of Reference vs Mutated Amino Acids HOPX Alphamissense", 
       x = "Reference Amino Acid", 
       y = "Mutated Amino Acid", 
       fill = "Pathogenic score") + 
  theme_minimal() +  # Minimal theme for cleaner look
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.background = element_rect(fill = 'black'), 
        panel.grid = element_blank()) +
  geom_point(data = merged[is.na(Germline_classification)==F], aes(shape=Germline_classification), colour='black', size =3.5) 