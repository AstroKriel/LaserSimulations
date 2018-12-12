# INITIALISE WORKSPACE
# ---------------
rm(list = ls())
graphics.off()

library(ggplot2)
library(gridExtra)

# ---------------
# Parameter: Load data
# ---------------
tau_r        = 20
var_file_num = 115
var_start    = 1500
var_end      = 2000


# ---------------
# LOAD AND PLOT DATA:
# ---------------
setwd(paste("~/Documents/University/Year3Sem2/FR_RESEARCH/SW-NK/test_SYS4_tau=", tau_r, sep = ""))

var_step  = 1
var_start = 50000
var_end   = 50500


# Plot Bifurcation Data
y_gamma   <- scan(file = "bif_gamma.txt")
y_gamma   <- y_gamma[y_gamma < 0.08]
y_extrema <- scan(file = "bif_extrema.txt")
y_extrema <- y_extrema[1:length(y_gamma)]
ggsave(filename = "Rplot_Bif_Diag.png", 
       ggplot(data.frame(x = y_gamma, 
                         y = y_extrema),
              aes(x = y_gamma,
                  y = y_extrema)) + geom_point(size = 0.1) + 
         labs(x = expression(eta), y = "Output Power", 
              title = bquote("Bifurcation Diagram: " ~ tau[R] == .(tau_r))) +
         theme(plot.title = element_text(hjust = 0.5)),
       width = 20, height = 12, dpi = 600, units = "cm", device='png')


# Plot Horizontal and Vertical Data
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
if (!all(y_hor == y_ver) && !file.exists(paste("Rplot_", var_file_num, "_hor_polar.png", sep = ""))) {
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


# Plot FFT data
y_FFT <- scan(paste("Data_",var_file_num,"_FFT.txt", sep = ""))
y_FFT <- y_FFT[c(seq(1, (length(y_FFT)/2), var_step))]
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
