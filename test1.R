
# =========================
# Local Strain Field Script
# =========================

# ---- Load Libraries ----
library(dplyr)
library(akima)      # For interpolation
library(pracma)     # For gradient calculation
library(ggplot2)

# ---- Load or Simulate Data ----
# Simulate some example data (you can replace this with your own CSV or dataframe)
set.seed(123)
n_points <- 100
time_points <- c(0, 1)

# Generate synthetic point positions
# tracks <- data.frame(
#   id = rep(1:n_points, each = 2),
#   time = rep(time_points, times = n_points),
#   x = rep(runif(n_points, 0, 10), each = 2),
#   y = rep(runif(n_points, 0, 10), each = 2)
# )


# Add small synthetic displacement for time = 1
# tracks <- tracks %>%
#   group_by(id) %>%
#   mutate(
#     x = ifelse(time == 1, x + 0.1 * x + rnorm(1, 0, 0.02), x),
#     y = ifelse(time == 1, y + 0.05 * y + rnorm(1, 0, 0.02), y)
#   ) %>%
#   ungroup()

tracks <- read_csv("sample1.csv")

# ---- Select Reference and Deformed Time Points ----
ref_time <- 0
def_time <- 1

ref <- tracks %>% filter(time == ref_time)
def <- tracks %>% filter(time == def_time)

# Merge positions by ID
merged <- inner_join(ref, def, by = "id", suffix = c("_ref", "_def"))

# Calculate displacements
merged <- merged %>%
  mutate(
    u = x_def - x_ref,
    v = y_def - y_ref
  )

# ---- Interpolate to Grid ----
x_seq <- seq(min(merged$x_ref), max(merged$x_ref), length.out = 100)
y_seq <- seq(min(merged$y_ref), max(merged$y_ref), length.out = 100)

u_interp <- with(merged, interp(x_ref, y_ref, u, xo = x_seq, yo = y_seq, linear = TRUE))
v_interp <- with(merged, interp(x_ref, y_ref, v, xo = x_seq, yo = y_seq, linear = TRUE))

# ---- Compute Gradients and Strains ----
dx <- diff(u_interp$x[1:2])
dy <- diff(u_interp$y[1:2])

# Compute gradients assuming unit spacing
grad_u <- gradient(u_interp$z)
grad_v <- gradient(v_interp$z)

# Scale by actual grid spacing
du_dx <- grad_u$X / dx
du_dy <- grad_u$Y / dy
dv_dx <- grad_v$X / dx
dv_dy <- grad_v$Y / dy

epsilon_xx <- du_dx
epsilon_yy <- dv_dy
epsilon_xy <- 0.5 * (du_dy + dv_dx)

# ---- Visualization Helper ----
plot_strain_field <- function(x, y, field, title = "Strain", palette = terrain.colors) {
  filled.contour(
    x, y, field,
    color.palette = palette,
    plot.title = title(main = title, xlab = "x", ylab = "y"),
    plot.axes = { axis(1); axis(2) }
  )
}

# ---- Plot Strain Fields ----
plot_strain_field(u_interp$x, u_interp$y, epsilon_xx, title = "Strain εxx")
plot_strain_field(u_interp$x, u_interp$y, epsilon_yy, title = "Strain εyy")
plot_strain_field(u_interp$x, u_interp$y, epsilon_xy, title = "Strain εxy")