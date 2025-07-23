
# =========================
# Local Strain Field Script
# =========================

# ---- Load Libraries ----
library(dplyr)
library(akima)      # For interpolation
library(pracma)     # For gradient calculation
library(ggplot2)

# ---- Load Data ----
# The example script originally generated synthetic data for two time points.
# To work with real tracking data containing many time points we read the CSV
# file provided with this repository.  It must contain four columns:
# `cell_id`, `frame`, `x` and `y`.
tracks <- read.csv("sample1.csv")

# Rename columns to generic names used later in the script
colnames(tracks) <- c("id", "time", "x", "y")

# Ensure ordering by time
tracks <- tracks %>% arrange(time, id)

# ---- Select Reference and Deformed Time Points ----
# The reference frame is typically the first (time = 0).  Strain fields will be
# computed for all later frames listed in `def_times`.
ref_time <- 0
def_times <- sort(unique(tracks$time))
def_times <- def_times[def_times > ref_time]

# ---- Strain Calculation Helper ----
compute_strain <- function(data, ref_time, def_time, grid_n = 100) {
  ref <- data %>% filter(time == ref_time)
  def <- data %>% filter(time == def_time)

  merged <- inner_join(ref, def, by = "id", suffix = c("_ref", "_def")) %>%
    mutate(
      u = x_def - x_ref,
      v = y_def - y_ref
    )

  x_seq <- seq(min(merged$x_ref), max(merged$x_ref), length.out = grid_n)
  y_seq <- seq(min(merged$y_ref), max(merged$y_ref), length.out = grid_n)

  u_interp <- with(merged, interp(x_ref, y_ref, u, xo = x_seq, yo = y_seq, linear = TRUE))
  v_interp <- with(merged, interp(x_ref, y_ref, v, xo = x_seq, yo = y_seq, linear = TRUE))

  dx <- diff(u_interp$x[1:2])
  dy <- diff(u_interp$y[1:2])

  grad_u <- gradient(u_interp$z)
  grad_v <- gradient(v_interp$z)

  du_dx <- grad_u$X / dx
  du_dy <- grad_u$Y / dy
  dv_dx <- grad_v$X / dx
  dv_dy <- grad_v$Y / dy

  list(
    x = u_interp$x,
    y = u_interp$y,
    epsilon_xx = du_dx,
    epsilon_yy = dv_dy,
    epsilon_xy = 0.5 * (du_dy + dv_dx)
  )
}

# ---- Visualization Helper ----
plot_strain_field <- function(x, y, field, title = "Strain", palette = terrain.colors) {
  filled.contour(
    x, y, field,
    color.palette = palette,
    plot.title = title(main = title, xlab = "x", ylab = "y"),
    plot.axes = { axis(1); axis(2) }
  )
}

# ---- Loop Over Time Points ----
for (t in def_times) {
  strain <- compute_strain(tracks, ref_time, t)

  plot_strain_field(strain$x, strain$y, strain$epsilon_xx,
                    title = paste0("Strain εxx (t=", t, ")"))
  plot_strain_field(strain$x, strain$y, strain$epsilon_yy,
                    title = paste0("Strain εyy (t=", t, ")"))
  plot_strain_field(strain$x, strain$y, strain$epsilon_xy,
                    title = paste0("Strain εxy (t=", t, ")"))
}
