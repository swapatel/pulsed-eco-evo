# Written by Abigail Adjei

#install and load exponential matrix package
install.packages("expm")
library(expm)


# Function for one persistence-boundary panel

make_panel <- function(alpha12,
                       alpha23,
                       lambda,
                       panel_label,
                       fixed_values = c(0.1, 0.3, 0.5), # delta_3s
                       n = 150,
                       tau = 20) {
  
  # Natural mortality
  mu <- 0.1
  

  # 3-stage ecological matrix
  
  A <- matrix(c(-(alpha12 + mu),0, lambda,
    alpha12,-(alpha23 + mu), 0,
    0, alpha23,-mu
  ),
  nrow = 3,
  byrow = TRUE)
  
  # Matrix exponential
  exp_A <- expm(tau * A)
  
  # Grid for delta_1 and delta_2
  grid <- seq(1/n, 1, length.out = n)
  
  # Same colors in every panel
  curve_cols <- c("red", "blue", "forestgreen")
  
  # Empty plot
  plot(
    c(1/n, 1),
    c(1/n, 1),
    type = "n",
    xaxs = "i",
    yaxs = "i",
    xlab = expression(delta[1]),
    ylab = expression(delta[2]),
    main = expression(r(D * e^{tau*A}) == 1),
    # Larger font sizes
    cex.axis = 1.25,
    cex.lab = 1.5,
    cex.main = 1.35
  )
  

  # Calculate the r = 1 contours while delta_1 and delta_2 vary
  for (k in seq_along(fixed_values)) {
    
    r_mat <- matrix(NA, nrow = n, ncol = n)
    
    for (i in 1:n) {
      for (j in 1:n) {
        
        delta1 <- grid[i]
        delta2 <- grid[j]
        delta3 <- fixed_values[k]
        
        # Pulse survival matrix
        D <- diag(c(delta1, delta2, delta3))
        
        # One-pulse map
        C <- D %*% exp_A
        
        # Spectral radius
        r_mat[i, j] <- max(abs(eigen(C)$values))
      }
    }
    
    # Draw threshold r = 1
    contour(
      x = grid,
      y = grid,
      z = r_mat,
      levels = 1,
      add = TRUE,
      drawlabels = FALSE,
      col = curve_cols[k],
      lwd = 3
    )
  }
  

  # Panel label
  mtext(
    panel_label,
    side = 3,
    line = 0.15,
    adj = 0,
    font = 2,
    cex = 1.35
  )
}


# =========================================================
# FOUR-PANEL FIGURE
# =========================================================

# Layout:
#
#       A       B       shared legend
#       C       D       shared legend
#
layout(
  matrix(
    c(1, 2, 5,
      3, 4, 5),
    nrow = 2,
    byrow = TRUE
  ),
  widths = c(1, 1, 0.42)
)

# Margins and square plotting regions
par(
  mar = c(4.5, 4.5, 3.5, 1.2),
  pty = "s"
)


# (A)
make_panel(
  alpha12 = 0.20,
  alpha23 = 0.45,
  lambda = 0.35,
  panel_label = "(A)"
)

# (B)
make_panel(
  alpha12 = 0.45,
  alpha23 = 0.20,
  lambda = 0.35,
  panel_label = "(B)"
)


# (C)
make_panel(
  alpha12 = 0.20,
  alpha23 = 0.45,
  lambda = 0.40,
  panel_label = "(C)"
)


# (D)
make_panel(
  alpha12 = 0.45,
  alpha23 = 0.20,
  lambda = 0.40,
  panel_label = "(D)"
)



# Shared Legend
par(mar = c(0, 0, 0, 0), pty = "m")

plot.new()

legend(
  "center",
  legend = c(
    expression(delta[3] == 0.1),
    expression(delta[3] == 0.3),
    expression(delta[3] == 0.5)
  ),
  col = c("red", "blue", "forestgreen"),
  lwd = 3,
  cex = 1.25,
  bty = "n"
)