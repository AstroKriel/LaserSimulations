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
sys_name  = "PCF"  # 'PROF' 'PRPCF', 'PCF'
var_param_set = 2       # parameter set in {1, 2}
var_tau_r     = 50

# Bifurcation Parameters
var_bif_min   = 0
var_bif_max   = 0.95


# ---------------
# LOAD AND PLOT DATA:
# ---------------
# Load Forward Data
setwd(paste("~/Documents/University/Year3Sem2/FR_RESEARCH/Code-NK/SYS_", sys_name,"_Param=", var_param_set,
            "_tau=", var_tau_r, "_F", sep = ""))

bif_fft_F <- rev(scan(file = "fft_display.txt"))
bif_eta_F <- rev(scan(file = "bif_eta.txt"))
bif_mask_F = (var_bif_min < bif_eta_F) & (bif_eta_F < var_bif_max)
bif_eta_F <- bif_eta_F[bif_mask_F]
bif_y_F <- rev(scan(file = "bif_extrema.txt"))
bif_y_F <- bif_y_F[bif_mask_F]

setwd("~/Documents/University/Year3Sem2/FR_RESEARCH/Code-NK/Plots")
ggsave(filename = paste("Bif_FFT_SYS_", sys_name,"_Param=", var_param_set, "_tau=", var_tau_r, "_F",".png", sep = ""),
       ggplot(data.frame(x = bif_eta_F,
                         y = bif_y_F),
              aes(x,
                  y)) + geom_point(size = 0.001) +
         labs(x = expression(eta), y = "Output Power",
              title = bquote(paste("System ", .(sys_name)," Bifurcation Diagram F: ",
                                   sep = "") ~ tau[R] == .(var_tau_r))) +
         theme(plot.title = element_text(hjust = 0.01)),
       width = 20, height = 12, dpi = 600, units = "cm", device='png')
