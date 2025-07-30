# ============================
# Engineering Strain Example
# ============================
#
# This script computes the engineering strain field between the first and
# last frame of `sample1.csv`.  The result is visualised as a colour map with
# particle positions (translation corrected) overlaid.

library(dplyr)
library(readr)
library(akima)
library(pracma)
library(ggplot2)
library(scales)

# ---- Load and prepare tracks ----
tracks <- read_csv("sample1.csv",
                   col_names = c("id", "frame", "x", "y"),
                   show_col_types = FALSE)

ref_frame    <- min(tracks$frame)
target_frame <- max(tracks$frame)

# baseline and target data
ref <- filter(tracks, frame == ref_frame)
def <- filter(tracks, frame == target_frame)

# match particles and remove rigid translation
m <- inner_join(ref, def, by = "id", suffix = c("_ref", "_def"))
if (nrow(m) == 0) stop("No common particles between frames")

tx <- median(m$x_def - m$x_ref, na.rm = TRUE)
ty <- median(m$y_def - m$y_ref, na.rm = TRUE)

m <- mutate(m,
            u = (x_def - tx) - x_ref,
            v = (y_def - ty) - y_ref)

# ---- Interpolate displacement to a grid ----
grid_n <- 150
x_seq  <- seq(min(m$x_ref), max(m$x_ref), length.out = grid_n)
y_seq  <- seq(min(m$y_ref), max(m$y_ref), length.out = grid_n)

ui <- with(m, interp(x_ref, y_ref, u, xo = x_seq, yo = y_seq, linear = TRUE))
vi <- with(m, interp(x_ref, y_ref, v, xo = x_seq, yo = y_seq, linear = TRUE))

# ---- Engineering strain ε_xx ----
dx <- diff(ui$x[1:2])
grad_u <- gradient(ui$z)
eps_xx <- grad_u$X / dx

strain_df <- expand.grid(x = ui$x, y = ui$y) %>%
  mutate(eps_xx = as.vector(eps_xx))

# ---- Translation-corrected particle positions for all frames ----
ref0 <- ref %>% select(id, x_ref = x, y_ref = y)

disp_vs_ref <- tracks %>%
  inner_join(ref0, by = "id") %>%
  mutate(dx = x - x_ref, dy = y - y_ref)

trans <- disp_vs_ref %>%
  group_by(frame) %>%
  summarise(tx = median(dx, na.rm = TRUE),
            ty = median(dy, na.rm = TRUE), .groups = "drop")

tracks_corr <- tracks %>%
  left_join(trans, by = "frame") %>%
  mutate(x_corr = x - tx, y_corr = y - ty) %>%
  filter(frame <= target_frame)

# ---- Plot: strain field + particle positions ----

p <- ggplot() +
  geom_raster(data = strain_df, aes(x, -y, fill = eps_xx), interpolate = TRUE) +
  scale_fill_viridis_c(name = expression(epsilon[xx]),
                       labels = percent_format(accuracy = 0.01)) +
  geom_point(data = tracks_corr, aes(x_corr, -y_corr, color = frame),
             size = 0.6, alpha = 0.6) +
  coord_equal() +
  labs(title = sprintf("Engineering strain between frames %d and %d",
                       ref_frame, target_frame),
       x = "x (corrected)", y = "-y (corrected)") +
  guides(color = guide_colorbar(title = "Frame")) +
  theme_minimal(base_size = 12)

print(p)
