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
var_sys_name  = "PROF"  # 'PROF' 'PRPCF', 'PCF'
var_param_set = 2       # parameter set in {1, 2}
var_tau_r     = 20

# Bifurcation Parameters
var_bif_max   = 0.11
var_bif_min   = 0

# ---------------
# LOAD AND PLOT DATA:
# ---------------
# Load Forward Data
setwd(paste("~/Documents/GitHub/NK_FR/FunctionalisedCode/Param", var_param_set,
            "_SYS_", var_sys_name,"_tau=", var_tau_r, "_F", sep = ""))

bif_eta_F <- rev(scan(file = "bif_eta.txt"))
bif_mask_F = (var_bif_min < bif_eta_F) & (bif_eta_F < var_bif_max)
bif_eta_F <- bif_eta_F[bif_mask_F]
bif_y_F <- rev(scan(file = "bif_extrema.txt"))
bif_y_F <- bif_y_F[bif_mask_F]

ggsave(filename = paste("Rplot_", var_sys_name,"_Bif_tau=", var_tau_r,"_param=", var_param_set,".png", sep = ""), 
       ggplot(data.frame(x = bif_eta_F,
                         y = bif_y_F),
              aes(x,
                  y)) + geom_point(size = 0.001) +
         labs(x = expression(eta), y = "Output Power",
              title = bquote(paste("System ", .(var_sys_name)," Bifurcation Diagram F: ",
                                   sep = "") ~ tau[R] == .(var_tau_r))) +
         theme(plot.title = element_text(hjust = 0.01)),
       width = 20, height = 12, dpi = 600, units = "cm", device='png')


# Load Backward Data
setwd(paste("~/Documents/GitHub/NK_FR/FunctionalisedCode/Param", var_param_set,
            "_SYS_", var_sys_name,"_tau=", var_tau_r, "_B", sep = ""))

bif_eta_B <- rev(scan(file = "bif_eta.txt"))
bif_mask_B = (var_bif_min < bif_eta_B) & (bif_eta_B < var_bif_max)
bif_eta_B <- bif_eta_B[bif_mask_B]
bif_y_B <- rev(scan(file = "bif_extrema.txt"))
bif_y_B <- bif_y_B[bif_mask_B]

ggsave(filename = paste("Rplot_", var_sys_name,"_Bif_tau=", var_tau_r,"_param=", var_param_set,".png", sep = ""), 
       ggplot(data.frame(x = bif_eta_B,
                         y = bif_y_B),
              aes(x,
                  y)) + geom_point(size = 0.001) +
         labs(x = expression(eta), y = "Output Power",
              title = bquote(paste("System ", .(var_sys_name)," Bifurcation Diagram B: ",
                                   sep = "") ~ tau[R] == .(var_tau_r))) +
         theme(plot.title = element_text(hjust = 0.01)),
       width = 20, height = 12, dpi = 600, units = "cm", device='png')


# Plot Bifurcation Comparison
setwd(paste("~/Documents/GitHub/NK_FR/FunctionalisedCode", sep = ""))

ggsave(filename = paste("Rplot_", var_sys_name,"_Bif_Comp_tau=", var_tau_r,"_param=", var_param_set,".png", sep = ""),
       ggplot(data.frame(x = bif_eta_F,
                         y = bif_y_F),
              aes(x,y)) + 
         geom_point(size = 0.001, aes(color = "Forward")) +
         labs(x = expression(eta), y = "Output Power",
              title = bquote(paste("System ", .(var_sys_name)," Bifurcation Diagram: Sweep Direction Comparison: ",
                                   sep = "") ~ tau[R] == .(var_tau_r)),
              colour = c(expression(paste(eta, " Sweep Direction:", sep = "")))) +
         theme(plot.title = element_text(hjust = 0.01),
               legend.position="bottom",
               legend.spacing.y = unit(0, "mm"), 
               panel.border = element_rect(colour = "black", fill=NA),
               legend.background = element_blank(),
               legend.box.background = element_rect(colour = "black")) + 
         
         geom_point(data = data.frame(x = bif_eta_B, y = bif_y_B),
                    mapping = aes(x, y, colour = "Backward"), size = 0.001),
       width = 20, height = 12, dpi = 600, units = "cm", device='png')
