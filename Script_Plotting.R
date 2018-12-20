# INITIALISE WORKSPACE
# ---------------
rm(list = ls())
graphics.off()

library(ggplot2)
library(gridExtra)

# ---------------
# Parameters
# ---------------
# System Parameters
var_file_num  = 39  # which file to analyse
var_sys       = 4   # which system in [1, 4] \ {3}
var_param_set = 1   # parameter set in {1, 2}
tau_r         = 50

# Bifurcation Plot Parameters
plot_bif      = TRUE
var_bif_max   = 0.1
var_bif_min   = 0
var_bif_max2  = 0.1
var_bif_min2  = 0

# Frequency Plot Parameters
plot_FFT      = FALSE
var_step_FFT  = 1

# Polarisation Time Plot Parameters
plot_time     = FALSE
var_start     = 50000
var_end       = 50500


# ---------------
# LOAD AND PLOT DATA:
# ---------------
setwd(paste("~/Documents/GitHub/NK_FR/Param", var_param_set,"_SYS", var_sys,"_tau=", tau_r, sep = ""))

# Plot Bifurcation Data
if (plot_bif) {
  y_eta   <- scan(file = "bif_gamma.txt")
  y_eta   <- y_eta[(var_bif_min < y_eta) & (y_eta < var_bif_max)]
  y_extrema <- scan(file = "bif_extrema.txt")
  y_extrema <- y_extrema[(var_bif_min < y_eta) & (y_eta < var_bif_max)]
  ggsave(filename = "Rplot_Bif_Diag.png", 
         ggplot(data.frame(x = y_eta, 
                           y = y_extrema),
                aes(x = y_eta,
                    y = y_extrema)) + geom_point(size = 0.1) + 
           labs(x = expression(eta), y = "Output Power", 
                title = bquote("Bifurcation Diagram: " ~ tau[R] == .(tau_r))) +
           theme(plot.title = element_text(hjust = 0.01)),
         width = 20, height = 12, dpi = 600, units = "cm", device='png')
}

# y_eta     <- scan(file = "bif_gamma.txt")
# y_eta     <- y_eta[(var_bif_min2 < y_eta) & (y_eta < var_bif_max2)]
# y_extrema <- scan(file = "bif_extrema.txt")
# y_extrema <- y_extrema[1:length(y_eta)]
# qplot(y_eta, y_extrema, xlim = c(0.075, 0.085), ylim = c(0.3, 1.25))

# Time Plot
if (plot_time) {
  # Vert polar data
  y_ver <- (scan(paste("Data_", var_file_num, "_ver_polar.txt", sep = "")))[var_start:var_end]
  x_ver <- (1:length(y_ver))
  ggsave(filename = paste("Rplot_", var_file_num, "_ver_polar.png", sep = ""), 
         qplot(x_ver,
               y_ver,
               geom = "line",
               main = "Vertical Polarisation",
               xlab = "",
               ylab = "") +
           theme(plot.title = element_text(hjust = 0.5)),
         width = 20, height = 12, dpi = 600, units = "cm", device='png')
  
  # Horiz polar data
  y_hor <- (scan(paste("Data_", var_file_num, "_hor_polar.txt", sep = "")))[var_start:var_end]
  if (!all(y_hor == y_ver)) {
    # if horiz and vert data are not equal
    x_hor <- (1:length(y_hor))
    ggsave(filename = paste("Rplot_", var_file_num, "_hor_polar.png", sep = ""), 
           qplot(x_hor,
                 y_hor,
                 geom = "line",
                 main = "Horizontal Polarisation",
                 xlab = "",
                 ylab = "") +
             theme(plot.title = element_text(hjust = 0.5)),
           width = 20, height = 12, dpi = 600, units = "cm", device='png')
  }
}


# FFT Plot
if (plot_FFT) {
  y_FFT <- scan(paste("Data_",var_file_num,"_FFT.txt", sep = ""))
  y_FFT <- y_FFT[c(seq(1, (length(y_FFT)/2), var_step_FFT))]
  x_FFT <- 1:length(y_FFT)
  ggsave(filename = paste("Rplot_", var_file_num, "_FFT.png", sep = ""),
         qplot(x_FFT,
               y_FFT,
               geom = "line",
               main = "Frequency Spectrum",
               xlab = "",
               ylab = "") +
           theme(plot.title = element_text(hjust = 0.5)),
         width = 20, height = 12, dpi = 600, units = "cm", device='png')
}


# ---------------
# Done
# ---------------
print("done")