
library(ggplot2)
library(ggpubr)
library(dplyr)


colors_ten_families <-   c("Acaulosporaceae" = "#9e0142",
                           "Ambisporaceae" = "#d53e4f",
                           "Archaeosporaceae" = "#f46d43",
                           "Claroideoglomeraceae" = "#fdae61",
                           "Diversisporaceae" = "#fee08b",
                           "Geosiphonaceae" = "#fbf898",
                           "Gigasporaceae" = "#daeb7f",
                           "Glomeraceae" = "#abdda4",
                           "Glomeromycotina" = "#66c2a5",
                           "Paraglomeraceae" = "#3288bd")

colors_nine_families <-   c("Acaulosporaceae" = "#9e0142",
                           "Ambisporaceae" = "#d53e4f",
                           "Archaeosporaceae" = "#f46d43",
                           "Claroideoglomeraceae" = "#fdae61",
                           "Diversisporaceae" = "#fee08b",
                           "Gigasporaceae" = "#daeb7f",
                           "Glomeraceae" = "#abdda4",
                           "Glomeromycotina" = "#66c2a5",
                           "Paraglomeraceae" = "#3288bd")

colors_eight_families <-   c("Acaulosporaceae" = "#9e0142",
                            "Ambisporaceae" = "#d53e4f",
                            "Archaeosporaceae" = "#f46d43",
                            "Claroideoglomeraceae" = "#fdae61",
                            "Gigasporaceae" = "#daeb7f",
                            "Glomeraceae" = "#abdda4",
                            "Glomeromycotina" = "#66c2a5",
                            "Paraglomeraceae" = "#3288bd")


make_barplot <- function (rbind_df, upper_limit, color_vector, title) {
  
  bar <- 
    ggplot(rbind_df, aes(x=factor(Host, level = host_names), y=rel_abun, fill=Family, alpha=`In Both Habitats`)) +
    geom_bar(stat="identity", position="stack") +
    theme(axis.text.x = element_text(angle=90, size=12, vjust=0.5, 
                                     hjust=1, colour="black", face="italic"), 
          legend.box.margin = margin(0,0,0,1, "cm"), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_blank(), 
          axis.line = element_line(colour="black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          legend.title = element_text(size=16),
          legend.text = element_text(size=12, colour="black"),
          axis.text.y = element_text(size=12, colour="black")) +
    scale_y_continuous(limits = c(0, upper_limit), expand = c(0,0)) +
    theme(plot.margin = unit(c(0.5,0,0,0.5), "cm")) +
    labs(x = " ", y = " ", fill="Family") +
    coord_cartesian(clip = "off") +
    #geom_text(data=bar_sum_df, aes(x=Host, y=abund, label=count), inherit.aes = F, vjust=-1) +
    ggtitle(paste(title)) +
    theme(plot.title = element_text(margin = margin(0,0,30,0),
                                    hjust = 0.5,
                                    size = 16,
                                    face = "bold")) +
    #ggtitle("\n ") +
    #theme(plot.title=element_text(margin = margin(0,0,30,0))) +
    scale_fill_manual(values = color_vector) +
    guides(fill = guide_legend(order = 1),
           alpha = guide_legend(order = 2)) +
    scale_alpha_discrete(range = c(0.45, 1), breaks = c("Yes", "No")) 
  
}


noncore_barplot <- function (rbind_df, upper_limit, color_vector, title) {
  
  bar <- 
    ggplot(rbind_df, aes(x=Host, y=rel_abun, fill=Family)) +
    geom_bar(stat="identity", position="stack") +
    theme(axis.text.x = element_text(angle=90, size=12, vjust=0.5, 
                                     hjust=1, colour="black", face="italic"), 
          legend.box.margin = margin(0,0,0,1, "cm"), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_blank(), 
          axis.line = element_line(colour="black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          legend.title = element_text(size=16),
          legend.text = element_text(size=12, colour="black"),
          axis.text.y = element_text(size=12, colour="black")) +
    scale_y_continuous(limits = c(0, upper_limit), expand = c(0,0)) +
    theme(plot.margin = unit(c(0.5,0,0,0.5), "cm")) +
    labs(x = " ", y = " ", fill="Family") +
    coord_cartesian(clip = "off") +
    #geom_text(data=bar_sum_df, aes(x=Host, y=abund, label=count), inherit.aes = F, vjust=-1) +
    ggtitle(paste(title)) +
    theme(plot.title = element_text(margin = margin(0,0,30,0),
                                    hjust = 0.5,
                                    size = 16,
                                    face = "bold")) +
    #ggtitle("\n ") +
    #theme(plot.title=element_text(margin = margin(0,0,30,0))) +
    scale_fill_manual(values = color_vector) 
  #guides(fill = guide_legend(order = 1),
  #       alpha = guide_legend(order = 2)) +
  #scale_alpha_discrete(range = c(1, 0.5), breaks = c("Yes", "No")) 
  
}



two_barplots <- function (ro_bar, ak_bar) {
  
  two_bars <- ggarrange(ro_bar, 
                        ak_bar, 
                        labels = c("A)", "B)"), 
                        vjust = 4.2, 
                        ncol = 2, 
                        common.legend = T, 
                        legend = "right") +
    theme(plot.margin = margin(0,1,0,0, "cm"))
  
  two_barplots <- annotate_figure(two_bars,
                                 left = text_grob("Relative Abundance (%)",
                                                  rot = 90,
                                                  size = 16),
                                 bottom = text_grob("Host", 
                                                    size = 16, 
                                                    hjust = 3.4))
  
  return(two_barplots)
  
}




##########################################################################

## keeping this just for reference
# how to make a 2x2 barplot

# first arrange side by side plots for shared and shared ASVs to get
# common legends for both (different families are present in each grouping)

# both_shared <- ggarrange(ro_shared_bar_noncore,
#                          ak_shared_bar_noncore, 
#                          labels=c("A)", "B)"), 
#                          hjust=-0.5, 
#                          vjust = 5.5, 
#                          common.legend = TRUE, 
#                          legend = "right")
# 
# both_shared_core <- ggarrange(ro_shared_bar_notitle,
#                               ak_shared_bar_notitle,
#                               labels=c("C)", "D)"),
#                               hjust=-0.5,
#                               vjust = 1,
#                               common.legend = TRUE,
#                               legend = "right")
# 
# 
# 
# # combine the pre-arranged side by side plots to make a 2x2 plot
# all_shared_barplots <- ggarrange(both_shared,
#                                  both_shared_core,
#                                  ncol=1,
#                                  nrow=2,
#                                  heights = c(1.15, 1)) +
#   theme(plot.margin = margin(0,0,0.5,0.5, "cm"))
# 
# all_shared_barplots <- annotate_figure(all_shared_barplots,
#                                        left = text_grob("Relative Abundance (%)",
#                                                         rot = 90,
#                                                         size = 16,
#                                                         hjust = 0.2),
#                                        bottom = text_grob("Host", size = 16, hjust = 2.7, vjust=-1))
# 
# ggsave("../figures/2x2_shared_barplots.png", width=11, height=12)


