# Clear workspace and console
rm(list = ls())
cat("\014") 

library(tidyverse)

################################################################################

# Discrete

################################################################################

plotAUC <- NULL
methods <- c("Cor_EN", "Cor_RF", "None_EN", "Lit_EN", "Lit_RF")
methodNames <- c("Correlation,\nElasticNet",
                 "Correlation,\nRandom Forest",
                 "ElasticNet",
                 "Literature,\nElasticNet",
                 "Literature,\nRandom Forest")
for (m in 1:length(methods)){
  if (!str_detect(methods[m], "Lit")){
    load(paste0("E:/Thesis/PublicationScripts/AUCdf_",methods[m],".RData"))
    AUCdf$Method <- methodNames[m]
    AUCdf$FactorName <- c("Low education",
                          "Physical\ninactivity",
                          "Unhealthy\ndiet",
                          "Depression",
                          "Type II\ndiabetes",
                          "Heart disease",
                          "Sex\n(male)")
    plotAUC <- rbind.data.frame(plotAUC, AUCdf)
  } else{
    load(paste0("E:/Thesis/PublicationScripts/AUCdf_",methods[m],".RData"))
    AUCdf$Method <- methodNames[m]
    AUCdf$FactorName <- c("Low education",
                          "Physical\ninactivity",
                          "Unhealthy\ndiet",
                          "Depression",
                          "Type II\ndiabetes",
                          "Heart disease")
    plotAUC <- rbind.data.frame(plotAUC, AUCdf) 
  }

}

p <- ggplot(plotAUC) +
  geom_bar(aes(x = Method, y = AUC, fill = FactorName, alpha = Method),
           stat = "identity", position = position_dodge(), color = "black") +
  facet_grid(cols = vars(FactorName)) +
  guides(fill = "none") +
  coord_cartesian(ylim = c(0.5,1)) +
  ylab("AUROC") +
  scale_alpha_manual(values = c(0.2, 0.4,0.6,0.8,1)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold",
                                  color = "white"))

g <- ggplot_gtable(ggplot_build(p))

strips <- which(grepl('strip-', g$layout$name))

pal <- RColorBrewer::brewer.pal(n = 7, "Dark2")


for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white"
}

plot(g)
ggsave(g, file = "Performance_discrete_EXTEND.png", width = 8, height = 5)

################################################################################

# Continuous

################################################################################

plotR2 <- NULL
methods <- c("Cor_EN", "Cor_RF", "None_EN", "Lit_EN", "Lit_RF")
methodNames <- c("Correlation,\nElasticNet",
                 "Correlation,\nRandom Forest",
                 "ElasticNet",
                 "Literature,\nElasticNet",
                 "Literature,\nRandom Forest")
for (m in 1:length(methods)){
  load(paste0("E:/Thesis/PublicationScripts/Data/MRS_Models/R2df_",methods[m],".RData"))
  R2df$Method <- methodNames[m]
  R2df$FactorName <- c("Systolic\nblood pressure",
                        "Total cholesterol")
  plotR2 <- rbind.data.frame(plotR2, R2df)
}

p <- ggplot(plotR2) +
  geom_bar(aes(x = Method, y = R2, fill = FactorName, alpha = Method),
           stat = "identity", position = position_dodge(), color = "black") +
  facet_grid(cols = vars(FactorName)) +
  guides(fill = "none") +
  ylab(expression(R^2)) +
  coord_cartesian(ylim = c(0.01,0.2)) +
  scale_alpha_manual(values = c(0.2, 0.4,0.6,0.8,1)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold",
                                  color = "white"))

g <- ggplot_gtable(ggplot_build(p))

strips <- which(grepl('strip-', g$layout$name))

pal <- RColorBrewer::brewer.pal(n = 7, "Set1")


for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white"
}

plot(g)

ggsave(g, file = "Performance_continuous1_EXTEND.png", width = 3, height = 5)
