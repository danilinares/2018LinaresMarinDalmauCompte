B <- 2000 # bootstrap samples
alpha <- .01
log_duration_min <- log10(0.005)
log_duration_max <- log10(0.3)

single_column_width <- 3.5
two_columns_width <- 7.2

size_point <- .7
size_point_cor <- 1.3
size_line <- .25
alpha_fill_level <- .25

label_duration <- "Duration (s)"
label_proportion <- "Proportion correct"
label_size <- "Size"
label_crt <- "Threshold on CRT (s)"
label_tablet <- "Threshold on tablet (s)"
label_platform <- "Platform"
label_small <- "Threshold for small (s)"
label_large <- "Threshold for large (s)"
label_crt_ss <- "Suppression index on CRT"
label_tablet_ss <- "Suppression index on tablet"

name_size <- c(`1` = "Small", `4` = "Large")

width_pix <- 800
height_pix <- 600
half_width_pix <- .5 * 800
half_height_pix <- .5 * 600
large_stimulus_size <- 111
small_stimulus_size <- 28


name_cond <- c("e" = "enceph", "c" = "control", "s" = "scz")

theme_set(theme_classic(base_size = 8) +
            theme(axis.line = element_line(size = size_line), 
                  axis.ticks = element_line(size = size_line), 
                  strip.background = element_rect(colour = "black", 
                                                  fill = "lightgrey",
                                                  linetype = 0)))
