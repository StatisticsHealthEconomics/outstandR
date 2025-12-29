library(ggplot2)
library(dplyr)
library(hexSticker)
library(showtext)
library(grid)

# --- Setup Fonts and Colors ---
# Load font (re-running just to be safe)
font_add_google("Montserrat", "montserrat")
showtext_auto()

c1 <- "#0072B2" # Blue
c2 <- "#D55E00" # Orange/Vermilion
custom_colors <- c("Group 1" = c1, "Group 2" = c2)


# --- 1. Data Generation (Same as before) ---
# Pie Chart
pie_data <- data.frame(group = factor(c("Group 1", "Group 2")), value = c(50, 50))
pie_plot <- ggplot(pie_data, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = custom_colors) +
  theme_void() + theme(legend.position = "none")
pie_grob <- ggplotGrob(pie_plot)

# Cloud Plot
set.seed(101)
cloud_data <- tibble(
  x = c(rnorm(70, mean = 3, sd = 1), rnorm(30, mean = 4.5, sd = 1)),
  y = c(rnorm(70, mean = 3, sd = 1), rnorm(30, mean = 4.5, sd = 1)),
  group = factor(c(rep("Group 1", 70), rep("Group 2", 30)))
)
outlier_x <- 7.5; outlier_y <- 7.5; outlier_size <- 1.5

p_main <- ggplot() +
  geom_point(data = cloud_data, aes(x = x, y = y, color = group),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors) +
  annotation_custom(grob = pie_grob,
                    xmin = outlier_x - outlier_size, xmax = outlier_x + outlier_size,
                    ymin = outlier_y - outlier_size, ymax = outlier_y + outlier_size) +
  xlim(-1, 11) + ylim(-1, 11) +
  theme_void() + theme_transparent() + theme(legend.position = "none")


# --- 2. Generate the Sticker Object (DO NOT SAVE YET) ---
# We set filename = NULL so it creates a ggplot object instead of a file.
sticker_object <- sticker(
  subplot = p_main,
  package = "outstandR",
  p_size = 22,
  # p_color = "#333333",
  p_color = "#D55E00",
  p_family = "montserrat",
  p_y = 0.6,
  s_x = 1, s_y = 1.15, s_width = 1.1, s_height = 1.1,
  h_fill = "#FFFFFF",
  h_color = c2,
  # We can use a nice thick border now because we will add a buffer
  h_size = 1.6
  # filename = NULL
)


# --- 3. The Buffer Hack ---
# Add a tiny margin around the plot object.
# This forces the canvas to be slightly larger than the hex border.
buffered_sticker <- sticker_object +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))

# --- 4. Save the Buffered Sticker ---
# Use ggsave explicitly. We set bg = "transparent" to keep the outside clear.
# The dimensions roughly approximate the standard hex ratio.
ggsave(filename = "outstandR_final_buffered.png",
       plot = buffered_sticker,
       width = 5.08, height = 5.86, units = "cm",
       bg = "transparent",
       dpi = 320)

showtext_auto(FALSE)
