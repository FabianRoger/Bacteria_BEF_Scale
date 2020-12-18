
# Title: Gamfeldt et al. Biodiversity and ecosystem functioning across gradients in spatial scale

# Project: customised plotting function

# create customised plotting theme
theme_meta <- function(base_size = 12, base_family = "") {
  theme(panel.background =  element_rect(fill = "white"), 
        panel.border = element_rect(fill="NA", color="black", size=0.35, linetype="solid"),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(),
        axis.ticks.length = unit(-0.16, "cm"),
        title = element_text(size = 7),
        axis.title.x = element_text(colour ="black", size = 9, face = "plain", margin=margin(5,0,0,0,"pt")),
        axis.title.y = element_text(colour = "black", size = 9, face = "plain", margin=margin(0,5,0,0,"pt")),
        axis.text.x = element_text(colour = "black", size=9, face = "plain",  margin=margin(10,0,0,0,"pt")),
        axis.text.y = element_text(colour ="black", size=9, face = "plain", margin=margin(0,10,0,0,"pt")),
        axis.ticks.x = element_line(colour = "black", size = 0.4),
        axis.ticks.y = element_line(colour = "black", size = 0.4),
        legend.title = element_text(size = 9, colour ="black", face = "plain"),
        legend.text = element_text(size = 9, colour ="black", face = "plain"))
}  

