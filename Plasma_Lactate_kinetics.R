"R Script for the visualisation of lactate kinetics in plasma samples of male and female control participants and professionals.
Output: Line plot"

library(tidyverse)
library(readxl)

data <- read_excel("/Users/majaewen/Documents/Uni/6_Semester/Daten/MSSA_TableS1-mmc1_02.01.2023.xlsx", sheet = 2)

data_long <- data %>%
  pivot_longer(cols = starts_with("Lactate"), 
               names_to = "Timepoint", 
               values_to = "Lactate") %>%
  mutate(Timepoint = case_when(
    Timepoint == "Lactate0" ~ "0 min",
    Timepoint == "Lactate7" ~ "7 min",
    Timepoint == "Lactate15" ~ "15 min",
    Timepoint == "LactatePeak" ~ "Peak",
    TRUE ~ Timepoint
  ))

data_long$Timepoint <- factor(data_long$Timepoint, levels = c("0 min", "7 min", "15 min", "Peak"))

farben <- c("CTR_M" = "blue", "PRO" = "red", "CTR_F" = "darkgrey")

# generate Line plot
p <- ggplot(data_long, aes(x = Timepoint, y = Lactate, group = Participant.ID, color = category)) +
  geom_line(alpha = 0.5) +       
  geom_point() +
  #stat_summary(aes(group = category), fun = mean, geom = "line", size = 1.5) + 
  scale_color_manual(values = farben) +
  labs(title = "Lactate kinetics", y = "Lactate (mmol/L)", x = " ") +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none") 

ggsave("Lactate kinetics.png", plot = p, width = 8, height = 6)