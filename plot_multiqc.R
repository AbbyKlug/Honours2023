library(ggplot2)
sample_names <- c("AD4", "AD5", "AD7", "AD9", "Old7", "Old8", "Old9","Old10", "Young1", "Young3", "Young5", "Young8")
map_perc <- c(52.4057, 49.3639, 52.9756, 49.3304, 53.2870, 57.3342, 56.8442, 43.5227, 68.2766, 65.4710, 53.3409, 46.1667)

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc)

# Bar plot

mapping_rate <- ggplot(data, aes(x = Sample, y = MappedPercentage)) +
  geom_bar(stat = "identity", fill = "hotpink") +
  geom_text(aes(label = MappedPercentage), angle = -90, vjust = 1, color = "black") +
  labs(
       x = "Sample", y = "Mapped Percentage (%)") +
  theme_minimal() +
  ylim(0, 100) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = -90, hjust = 1))  # Rotate x-axis labels

mapping_rate <- edit_plots(mapping_rate)
mapping_rate

ggplot(data, aes(x = Sample, y = MappedPercentage)) +
  geom_bar(stat = "identity", fill = "hotpink") +
  geom_text(aes(label = MappedPercentage), angle = -90, vjust = 1, color = "black") +
  labs(
    x = "Sample", y = "Mapped Percentage (%)") +
  theme_minimal() +
  ylim(0, 100) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = -90, hjust = 1))  # Rotate x-axis labels

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
mapping_rate <- edit_plots(mapping_rate)
mapping_rate


ggplot(data, aes(x = Sample, y = MappedPercentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = MappedPercentage), vjust = -0.5, color = "black") +
  labs(title = "Percentage of Mapped Reads by Sample",
       x = "Sample", y = "Mapped Percentage (%)") +
  theme_minimal()
