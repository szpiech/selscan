library(ggplot2)
library(reshape2)

# --- Read the data ---
ehh_data <- read.table("outfile.ehh.1943.out", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Rename derEHH to derived for plotting
colnames(ehh_data)[colnames(ehh_data) == "derEHH"] <- "derived"

# Choose x-axis: pdist or gdist
x_axis <- "pdist"  # Change to "gdist" if needed

ehh_data$pdist <- ehh_data$pdist / 1e6

# Reshape data to long format
ehh_long <- melt(ehh_data, id.vars = x_axis, 
                 measure.vars = c("derived", "ancEHH", "EHH"),
                 variable.name = "Type", value.name = "EHH_Value")

# Define custom colors
ehh_colors <- c("derived" = "red", "ancEHH" = "blue", "EHH" = "gray50")

# Plot
ggplot(ehh_long, aes_string(x = x_axis, y = "EHH_Value", color = "Type")) +
  geom_line(size = 1) +
  scale_color_manual(values = ehh_colors) +
  labs(x = ifelse(x_axis == "pdist", "Physical Distance (bp)", "Genetic Distance (cM)"),
       y = "EHH",
       title = "EHH Decay Plot") +
  theme_minimal() +
  theme(text = element_text(size = 14)) 