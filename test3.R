# ===========================================================
# ε_xx field and average ε_xx per frame, using known schedule
# Input: tracks(id, frame, x, y)
# Output:
#   - eps_xx_maps[[k]]: matrix of ε_xx(x,y) for frame k
#   - eps_summary: data.frame with global and field-avg ε_xx per frame
# ===========================================================

library(dplyr)
library(MASS)     # robust regression (rlm)
library(akima)    # interp to grid
library(pracma)   # gradient (assumes unit spacing)

tracks <- read_csv("sample1.csv")
# tracks <- filter(tracks, frame %in% seq(0, max(tracks$frame), 4))
# Rename columns to generic names used later in the script
colnames(tracks) <- c("id", "frame", "x", "y")

# Ensure ordering by time
# tracks <- tracks %>% arrange(frame, id)

# ------------------ Parameters ------------------
zero_frames <- 40          # frames 0..9: no applied strain
d_eps_per_frame <- 0.005/4   # +0.5% engineering strain per frame
max_total_eps <- 0.50      # up to 50% total
grid_n <- 120              # grid resolution for maps
use_robust <- TRUE         # robust fit (Huber) vs least squares
frames <- sort(unique(tracks$frame))

# -------------- Helpers (schedule & fit) --------------
lambda_x_inc <- function(k, zero_frames = 10, d_eps = 0.005, max_eps = 0.50) {
  # multiplicative axial stretch from k-1 → k
  if (k <= zero_frames) return(1.0)
  total_eps_k  <- min(max(0, (k - zero_frames) * d_eps), max_eps)
  total_eps_km1 <- min(max(0, (k - 1 - zero_frames) * d_eps), max_eps)
  if (total_eps_k == total_eps_km1) return(1.0)
  1.0 + d_eps
}

# Constrained affine: fix a11 = lambda_x, solve for a12, a21, a22, c1, c2
affine_fit_a11_fixed <- function(X, Y, a11, robust = TRUE) {
  x <- X[,1]; y <- X[,2]
  xp <- Y[,1]; yp <- Y[,2]
  
  A1 <- cbind(y, 1)          # for x': xp - a11*x = a12*y + c1
  b1 <- xp - a11 * x
  A2 <- cbind(x, y, 1)       # for y': yp = a21*x + a22*y + c2
  b2 <- yp
  
  if (robust) {
    p1 <- coef(MASS::rlm(A1, b1, psi = MASS::psi.huber, maxit = 200))
    p2 <- coef(MASS::rlm(A2, b2, psi = MASS::psi.huber, maxit = 200))
  } else {
    p1 <- qr.solve(A1, b1); p2 <- qr.solve(A2, b2)
  }
  
  a12 <- p1[1]; c1 <- p1[2]
  a21 <- p2[1]; a22 <- p2[2]; c2 <- p2[3]
  F <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
  cvec <- c(c1, c2)
  list(F = F, c = cvec)
}

# -------------- Main loop over frames --------------
eps_xx_maps <- list()  # ε_xx field per frame (matrix on a regular grid)
eps_summary <- list()  # per-frame summary numbers
Lambda_cum <- 1.0      # cumulative axial stretch from start (multiplicative)

for (k in frames[-1]) {
  k_prev <- k - 1
  ref <- tracks %>% filter(frame == k_prev)
  def <- tracks %>% filter(frame == k)
  m <- inner_join(ref, def, by = "id", suffix = c("_ref", "_def"))
  if (nrow(m) < 6) next
  
  # Re-center around median to numerically stabilize the fit
  x0 <- median(m$x_ref); y0 <- median(m$y_ref)
  X <- cbind(m$x_ref - x0, m$y_ref - y0)
  Y <- cbind(m$x_def - x0, m$y_def - y0)
  
  lamx <- lambda_x_inc(k, zero_frames, d_eps_per_frame, max_total_eps)
  fit <- affine_fit_a11_fixed(X, Y, a11 = lamx, robust = use_robust)
  
  # Predicted bulk motion and residuals
  Y_pred <- (X %*% t(fit$F)) + matrix(fit$c, nrow(X), 2, byrow = TRUE)
  U_res <- Y - Y_pred
  
  # ---- Interpolate residual u,v to a grid ----
  res_df <- data.frame(x = m$x_ref, y = m$y_ref,
                       u_res = U_res[,1], v_res = U_res[,2])
  
  x_seq <- seq(min(res_df$x), max(res_df$x), length.out = grid_n)
  y_seq <- seq(min(res_df$y), max(res_df$y), length.out = grid_n)
  
  ui <- with(res_df, akima::interp(x, y, u_res, xo = x_seq, yo = y_seq, linear = TRUE))
  vi <- with(res_df, akima::interp(x, y, v_res, xo = x_seq, yo = y_seq, linear = TRUE))
  
  dx <- diff(ui$x[1:2]); dy <- diff(ui$y[1:2])
  gu <- pracma::gradient(ui$z); gv <- pracma::gradient(vi$z)
  
  du_dx <- gu$X / dx
  # du_dy <- gu$Y / dy   # only needed if you also want shear
  # dv_dx <- gv$X / dx
  # dv_dy <- gv$Y / dy
  
  # ---- ε_xx field (incremental) ----
  eps_xx_residual <- du_dx                     # residual part from local non-uniformity
  eps_xx_global   <- fit$F[1,1] - 1            # equals +0.5% after frame 10 by construction
  eps_xx_total    <- eps_xx_global + eps_xx_residual  # small-increment composition
  
  eps_xx_maps[[as.character(k)]] <- list(x = ui$x, y = ui$y, z = eps_xx_total)
  
  # ---- averages ----
  avg_field <- mean(eps_xx_total, na.rm = TRUE) # area-average of the map
  avg_global <- eps_xx_global                   # cleanest estimate of average ε_xx
  
  # ---- cumulative axial strain (from start) ----
  Lambda_cum <- Lambda_cum * fit$F[1,1]            # multiplicative accumulation
  eng_strain_total <- Lambda_cum - 1               # engineering strain to date
  true_strain_total <- log(Lambda_cum)             # "true"/log strain
  
  eps_summary[[as.character(k)]] <- data.frame(
    frame = k,
    eps_xx_inc_global = eps_xx_global,
    eps_xx_inc_field_mean = avg_field,
    Lambda_cum = Lambda_cum,
    eps_xx_total_engineering = eng_strain_total,
    eps_xx_total_true = true_strain_total
  )
}

eps_summary <- dplyr::bind_rows(eps_summary)

# ----------------- Example: plot ε_xx field for the last frame -----------------
if (length(eps_xx_maps) > 0) {
  lastk <- as.integer(names(eps_xx_maps))[length(eps_xx_maps)]
  m <- eps_xx_maps[[as.character(lastk)]]
  filled.contour(m$x, m$y, m$z,
                 color.palette = terrain.colors,
                 plot.title = title(main = sprintf("ε_xx field (frame %d)", lastk),
                                    xlab = "x", ylab = "y"),
                 key.title = title("asdffffasdf"), omd = seq(10,4),
                 plot.axes = { axis(1); axis(2) }
                 )
}

plot_strain_field <- function(x, y, field, title = expression(epsilon[xx]),
                              nbreaks = 8, show_percent = TRUE,
                              right_margin = 7, cex_key = 0.85) {
  zrange <- range(field, na.rm = TRUE)
  levs <- pretty(zrange, n = nbreaks)
  lab <- if (show_percent) sprintf("%.2f%%", 100*levs) else sprintf("%.3f", levs)
  
  filled.contour(
    x, y, field,
    levels = levs,
    col = hcl.colors(length(levs) - 1, "YlGnBu", rev = TRUE),
    mar = c(4.5, 4.5, 2, right_margin),        # <- more space for the key labels
    plot.title = title(main = title, xlab = "x", ylab = "y"),
    plot.axes  = { axis(1); axis(2) },
    key.axes   = { axis(4, at = levs, labels = lab, las = 1, cex.axis = cex_key) }
  )
}
plot_strain_field(ui$x, ui$y, eps_xx_total, title = "ε_xx (total)", show_percent = F)
plot_strain_field(xg, yg, eps_xx_total, show_percent = TRUE)

library(ggplot2)

df <- expand.grid(x = ui$x, y = ui$y)

df$eps <- as.vector(eps_xx_total)
max(df$eps, na.rm = T)

ggplot(df, aes(x, y, fill = eps)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(name = expression(epsilon[xx]), labels = scales::percent_format(accuracy = 0.01)) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")


# ----------------- Inspect averages -----------------
print(head(eps_summary))
eps_summary
filter(eps_summary, frame %in% seq(0, max(tracks$frame), 8))


# disp plot 1 ----
library(dplyr)
library(tidyr)
library(ggplot2)

# tracks must have: id, frame, x, y
ref_frame <- min(tracks$frame)   # or set to 0 if you know that’s your baseline

# 1) Displacement relative to reference frame (per particle)
ref0 <- tracks %>% filter(frame == ref_frame) %>% dplyr::select(id, x0 = x, y0 = y)
d <- tracks %>%
  inner_join(ref0, by = "id") %>%
  mutate(ux = x - x0, vy = y - y0)

# 2) Robust global translation per frame (median across all particles)
trans <- d %>%
  group_by(frame) %>%
  summarise(tx = median(ux, na.rm = TRUE),
            ty = median(vy, na.rm = TRUE),
            .groups = "drop")

# 3) Translation-corrected displacements
d_corr <- d %>%
  left_join(trans, by = "frame") %>%
  mutate(u_res = ux - tx,            # x-displacement after translation removal
         v_res = vy - ty,            # y-displacement after translation removal
         disp_res = sqrt(u_res^2 + v_res^2))

# 4) Long format for facetting different components
d_long <- d_corr %>%
  dplyr::select(id, frame, u_res, v_res, disp_res) %>%
  pivot_longer(c(u_res, v_res, disp_res),
               names_to = "component", values_to = "value")

# 5) Plot: per-particle trajectories (faint) + median overlay
ggplot(d_long, aes(frame, value, group = interaction(id, component))) +
  geom_line(alpha = 0.15) +
  stat_summary(aes(group = component), fun = median, geom = "line", size = 1.1) +
  facet_wrap(~ component, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(u_res = "u (x, translation-corrected)",
                                      v_res = "v (y, translation-corrected)",
                                      disp_res = "|u| (magnitude, corrected)"))) +
  geom_vline(xintercept = 10, linetype = 2) +  # strain starts at frame 10
  theme_minimal(base_size = 12) +
  labs(title = "Particle displacement over time (translation-corrected)",
       x = "Frame", y = "Displacement (same units as x,y)")


# disp plot 2 ----
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Per-frame increments for each particle
inc <- tracks %>%
  arrange(id, frame) %>%
  group_by(id) %>%
  mutate(dx = x - dplyr::lag(x),
         dy = y - dplyr::lag(y)) %>%
  ungroup() %>%
  filter(!is.na(dx))

# 2) Per-frame global translation increment (robust)
trans_inc <- inc %>%
  group_by(frame) %>%
  summarise(tx = median(dx, na.rm = TRUE),
            ty = median(dy, na.rm = TRUE),
            .groups = "drop")

# 3) Residual increments and cumulative, per particle
inc_corr <- inc %>%
  left_join(trans_inc, by = "frame") %>%
  mutate(dx_res = dx - tx,
         dy_res = dy - ty) %>%
  group_by(id) %>%
  mutate(u_res = cumsum(coalesce(dx_res, 0)),
         v_res = cumsum(coalesce(dy_res, 0)),
         disp_res = sqrt(u_res^2 + v_res^2)) %>%
  ungroup()

# 4) Long format and plot
d_long <- inc_corr %>%
  dplyr::select(id, frame, u_res, v_res, disp_res) %>%
  pivot_longer(c(u_res, v_res, disp_res),
               names_to = "component", values_to = "value")

ggplot(d_long, aes(frame, value, group = interaction(id, component))) +
  geom_line(alpha = 0.15) +
  stat_summary(aes(group = component), fun = median, geom = "line", size = 1.1) +
  facet_wrap(~ component, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(u_res = "u (x, translation-corrected, cumulative)",
                                      v_res = "v (y, translation-corrected, cumulative)",
                                      disp_res = "|u| (magnitude, cumulative)"))) +
  geom_vline(xintercept = 10, linetype = 2) +
  theme_minimal(base_size = 12) +
  labs(title = "Particle displacement over time (translation-corrected, incremental accumulation)",
       x = "Frame", y = "Displacement (same units as x,y)")


# plot 3 ----
library(dplyr)
library(ggplot2)

# tracks must have: id, frame, x, y
ref_frame <- min(tracks$frame)   # or set explicitly, e.g. ref_frame <- 0

# --- 1) Compute per-frame translation relative to the reference frame ---
ref0 <- tracks %>% filter(frame == ref_frame) %>% dplyr::select(id, x_ref = x, y_ref = y)

# Join each frame to the reference by id, get displacements
disp_vs_ref <- tracks %>%
  inner_join(ref0, by = "id") %>%
  mutate(dx = x - x_ref,
         dy = y - y_ref)

# Robust translation estimate per frame (median across particles)
trans_by_frame <- disp_vs_ref %>%
  group_by(frame) %>%
  summarise(tx = median(dx, na.rm = TRUE),
            ty = median(dy, na.rm = TRUE),
            .groups = "drop") %>%
  # Ensure the reference frame has zero translation
  mutate(tx = ifelse(frame == ref_frame, 0, tx),
         ty = ifelse(frame == ref_frame, 0, ty))

# --- 2) Apply translation correction to all points ---
tracks_corr <- tracks %>%
  left_join(trans_by_frame, by = "frame") %>%
  mutate(tx = coalesce(tx, 0),   # in case a frame had no matches (rare)
         ty = coalesce(ty, 0),
         x_corr = x - tx,
         y_corr = y - ty)

# --- 3) Plot: x vs -y with color = frame (all frames overlaid) ---
ggplot(tracks_corr, aes(x_corr, -y_corr, color = frame)) +
  geom_point(alpha = 0.7, size = 0.9) +
  coord_equal() +
  labs(title = "Particle positions (translation-corrected)",
       x = "x (translation-corrected)",
       y = "-y (translation-corrected)") +
  guides(color = guide_colorbar(title = "Frame")) +
  theme_minimal(base_size = 12)


# plot 4 ----
library(dplyr)
library(ggplot2)
library(scales)   # for percent_format
# If you need viridis colors (nice for fields): install.packages("viridis"); library(viridis)

# ---- Inputs you control ----
ref_frame <- min(tracks$frame)   # or set explicitly, e.g. ref_frame <- 0
max_frame <- 120                  # <-- plot points up to this frame (inclusive)

# ---- 1) Translation correction (relative to reference frame) ----
ref0 <- tracks %>%
  filter(frame == ref_frame) %>%
  select(id, x_ref = x, y_ref = y)

# Displacements vs reference, per point and frame
disp_vs_ref <- tracks %>%
  inner_join(ref0, by = "id") %>%
  mutate(dx = x - x_ref,
         dy = y - y_ref)

# Robust per-frame translation (median across all points)
trans_by_frame <- disp_vs_ref %>%
  group_by(frame) %>%
  summarise(tx = median(dx, na.rm = TRUE),
            ty = median(dy, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(tx = ifelse(frame == ref_frame, 0, tx),
         ty = ifelse(frame == ref_frame, 0, ty))

# Apply translation correction to all points
tracks_corr <- tracks %>%
  left_join(trans_by_frame, by = "frame") %>%
  mutate(tx = coalesce(tx, 0),
         ty = coalesce(ty, 0),
         x_corr = x - tx,
         y_corr = y - ty)

# ---- 2a) Prepare STRAIN background if you have ui$x, ui$y, eps_xx_total ----
# (This is the ε_xx map for frame = max_frame, on the same physical coordinates.)
# Align the grid to the same translation correction used for that frame:
tx_k <- trans_by_frame$tx[trans_by_frame$frame == max_frame]
ty_k <- trans_by_frame$ty[trans_by_frame$frame == max_frame]

# Example if your interpolation is named `ui` and the matrix is `eps_xx_total`:
grid_df <- expand.grid(x = ui$x, y = ui$y)
grid_df$eps_xx <- as.vector(eps_xx_total)

# ---- 2b) ...or if you stored it as eps_xx_maps[[k]] with fields x,y,z ----
# (Uncomment this block and comment out the 2a block above if you're using maps)
# m <- eps_xx_maps[[as.character(max_frame)]]
# grid_df <- expand.grid(x = m$x, y = m$y)
# grid_df$eps_xx <- as.vector(m$z)

# After creating grid_df with columns x, y, eps_xx, translate it:
grid_df <- grid_df %>%
  mutate(x_corr = x - tx_k,
         y_corr = y - ty_k)

# ---- 3) One plot: ε_xx background + translation-corrected particle positions ----
p <- ggplot() +
  # Background field (flip Y by plotting -y to match your previous convention)
  geom_raster(
    data = grid_df,
    aes(x_corr, -y_corr, fill = eps_xx),
    interpolate = TRUE,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(
    name = expression(epsilon[xx]),
    labels = percent_format(accuracy = 0.01)
  ) +
  # Points from all frames up to max_frame
  geom_point(
    data = tracks_corr %>% filter(frame <= max_frame),
    aes(x_corr, -y_corr, color = frame),
    size = 0.7,
    alpha = 0.6
  ) +
  coord_equal() +
  labs(
    title = sprintf("ε[x][x] background + particle positions up to frame %d (translation-corrected)", max_frame),
    x = "x (translation-corrected)",
    y = "-y (translation-corrected)"
  ) +
  guides(
    color = guide_colorbar(title = "Frame"),
    fill  = guide_colorbar(title = expression(epsilon[xx]))
  ) +
  theme_minimal(base_size = 12)

print(p)



# full strain -----
library(dplyr)
library(akima)
library(pracma)
library(MASS)

# --- parameters ---
frame0 <- min(tracks$frame)   # baseline frame
k_bg   <- max(tracks$frame)   # frame to plot (set to your target)
k_bg <- 120 +80

# --- join baseline and target frame by id ---
ref <- tracks %>% filter(frame == frame0)
def <- tracks %>% filter(frame == k_bg)
m   <- inner_join(ref, def, by = "id", suffix = c("_ref", "_def"))

# --- optional: remove GLOBAL TRANSLATION only (robust median) ---
tx <- median(m$x_def - m$x_ref, na.rm = TRUE)
ty <- median(m$y_def - m$y_ref, na.rm = TRUE)

# translation-corrected displacement (baseline → frame k_bg)
m <- m %>%
  mutate(u = (x_def - tx) - x_ref,
         v = (y_def - ty) - y_ref)

# --- sanity check: expected vs observed global stretch ---
expected_eps <- min(max(0, k_bg - 40) * 0.005/4, 0.50)
message(sprintf("Expected ε_xx at frame %d ≈ %.1f%%", k_bg, 100*expected_eps))

# affine fit (unconstrained) to get global F between baseline and k_bg
X <- as.matrix(m[, c("x_ref", "y_ref")])
Y <- as.matrix(cbind(m$x_def - tx, m$y_def - ty))

affine_fit <- function(X, Y) {
  n <- nrow(X)
  M <- matrix(0, n*2, 6); b <- as.numeric(t(Y))
  M[seq(1, 2*n, 2), 1:2] <- X; M[seq(1, 2*n, 2), 3] <- 1
  M[seq(2, 2*n, 2), 4:5] <- X; M[seq(2, 2*n, 2), 6] <- 1
  p <- qr.solve(M, b)
  F <- matrix(c(p[1], p[2], p[4], p[5]), 2, 2, byrow = TRUE)
  c <- c(p[3], p[6]); list(F=F, c=c)
}
fit <- affine_fit(X, Y)
message("Observed global F between baseline and frame k_bg:\n"); print(fit$F)
message(sprintf("Observed engineering ε_xx (global) ≈ %.2f%%", 100*(fit$F[1,1]-1)))

# --- grid and gradients (TOTAL, baseline→k_bg) ---
x_seq <- seq(min(m$x_ref), max(m$x_ref), length.out = 150)
y_seq <- seq(min(m$y_ref), max(m$y_ref), length.out = 150)

ui <- with(m, akima::interp(x_ref, y_ref, u, xo = x_seq, yo = y_seq, linear = TRUE))
vi <- with(m, akima::interp(x_ref, y_ref, v, xo = x_seq, yo = y_seq, linear = TRUE))

dx <- diff(ui$x[1:2]); dy <- diff(ui$y[1:2])
gu <- pracma::gradient(ui$z); gv <- pracma::gradient(vi$z)
du_dx <- gu$X / dx; du_dy <- gu$Y / dy
dv_dx <- gv$X / dx; dv_dy <- gv$Y / dy

# Small-strain field (engineering) between baseline and k_bg:
eps_xx_small <- du_dx

# Finite strain (Green–Lagrange) if you prefer for large strains:
# F_field = I + Grad u  (Grad u wrt baseline coordinates)
F11 <- 1 + du_dx; F12 <- du_dy
F21 <- dv_dx;     F22 <- 1 + dv_dy
# E = 0.5*(F^T F - I);  E_xx:
E_xx <- 0.5 * (F11^2 + F21^2 - 1)

# Choose which to plot:
eps_to_plot <- eps_xx_small  # or E_xx for finite-strain measure

# --- ggplot overlay (field background + points up to k_bg) ---
library(ggplot2)
library(scales)

# prepare field data frame (also translation-correct grid to match your point plot convention)
grid_df <- expand.grid(x = ui$x, y = ui$y)
grid_df$eps <- as.vector(eps_to_plot)
grid_df <- grid_df %>% mutate(x_corr = x, y_corr = y)  # translation already applied to points below

# translation-corrected coordinates for points up to k_bg
tracks_corr <- tracks %>%
  filter(frame <= k_bg) %>%
  left_join(ref %>% select(id, x_ref=x, y_ref=y), by="id") %>%
  mutate(x_corr = x - tx, y_corr = y - ty)

ggplot() +
  geom_raster(data = grid_df, aes(x_corr, -y_corr, fill = eps), interpolate = TRUE, alpha = 0.9) +
  scale_fill_viridis_c(name = expression(epsilon[xx]), labels = percent_format(accuracy = 0.01)) +
  geom_point(data = tracks_corr, aes(x_corr, -y_corr, color = frame), size = 0.6, alpha = 0.6) +
  coord_equal() +
  labs(title = sprintf("Total ε[x][x] (baseline → frame %d) with translation-corrected points", k_bg),
       x = "x (corr.)", y = "-y (corr.)") +
  guides(color = guide_colorbar(title = "Frame")) +
  theme_minimal(base_size = 12)
