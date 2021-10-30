install.packages("BiocManager")
BiocManager::install(c("multtest","phyloseq","rhdf5","ggplot2","colorspace","stringi"))


library(ggplot2); library(phyloseq);library(vegan);library(multcomp);library(dplyr);library(grid);library(scales)
require(gridExtra); library(emmeans);library(multcompView); library(ggpubr); library(Rmisc); library(RColorBrewer)
library(purrr); library(RVAideMemoire); library(tibble); library(reshape2)

###################################
####Import files and pre-process###
###################################
#Import files generated in QIIME


unrarefied_physeq <- readRDS("data/unrarefied_physeq.rds")

otu_table(unrarefied_physeq) <- t(otu_table(unrarefied_physeq))

#kk_vt_rel_abund <- transform_sample_counts(unrarefied_physeq,function(x)x/sum(x))


# test <- tax_glom(unrarefied_physeq, "Genus")
# test_rel_abund <- transform_sample_counts(test, function(x)x/sum(x))
# test_melt <- psmelt(test_rel_abund)

#Melt phyloseq object to make a dataframe for ggplot and bipartite
#kk_vt.melt<-psmelt(kk_vt_rel_abund)
#names(kk_vt.melt)[1] <- "ASV"

kk_melt_counts <- psmelt(unrarefied_physeq)
names(kk_melt_counts)[1] <- "ASV"


plot_names <- sort(unique(kk_melt_counts$Plot))
host_names <- c("A. koa",
                "C. trigynum",
                "Grass",
                "M. polymorpha",
                "M. lessertiana",
                "R. hawaiensis",
                "Soil")


#######################
###Bipartite Network###
#######################
#load packages
library(sna);library(bipartite);library(Hmisc);library(ggpubr); library(Rmisc); library(emmeans);library(multcompView); 
library(multcomp); require(gridExtra); library(tidyverse)


#Convert a datatable into a network matrix for further use in bipartite.
#Create output as "list" type to calculate individual metrics and generate individual graphs for each plot.
#Prefereabble for different webs do not have to include all species names, i.e. they can be of different dimensions 
#(ragged). As such they are better suited for webs with non-comparable species sets (differing community rich/comps).
#Generates a weighted matrix for each plot


################################ Using Counts #####################################


#Convert a datatable into a network matrix for further use in bipartite.
#Create output as "list" type to calculate individual metrics and generate individual graphs for each plot.
#Prefereabble for different webs do not have to include all species names, i.e. they can be of different dimensions 
#(ragged). As such they are better suited for webs with non-comparable species sets (differing community rich/comps).
#Generates a weighted matrix for each plot

#kk_haka_bipartite <- frame2webs(kk_vt.melt,varnames=c("ASV","NewSampleType","Plot","Abundance"),type.out="list")

haka_counts_bipartite <- frame2webs(kk_melt_counts,
                                    varnames=c("ASV",
                                               "NewSampleType",
                                               "Plot",
                                               "Abundance"),
                                    type.out="list")

#haka_relabun <- haka_counts_bipartite

# remove soil
for (a_plot in plot_names) {
  
  c <- as.data.frame(haka_counts_bipartite[[a_plot]])
  
  #c <- c %>% dplyr::mutate(across(.cols = everything(), ~ .x/sum(.x)))
  
  c <- c[,1:6]
  
  c <- c[rowSums(c)>0,]
  
  c <- as.matrix(c)
  
  haka_counts_bipartite[[a_plot]] <- c
  
}


#kk_haka_composite_bipartite<-frame2webs(kk_vt.melt,varnames=c("ASV","NewSampleType","HabitatType","Abundance"),
#type.out="list")



##################################
###   Network specialization   ###
##################################
ct_network_specialization <- lapply(haka_counts_bipartite, networklevel, index= 'H2')

ct_H2_dat <- data.frame(unlist(ct_network_specialization))
row.names(ct_H2_dat) <- plot_names
ct_H2_dat$plot <- plot_names
ct_H2_dat$HabitatType <- rep(c("Restored","Remnant"),each=6)
names(ct_H2_dat)[1] <- "networkspecialization"


ct_H2_plot <- ggplot(ct_H2_dat,aes(x=HabitatType,y=networkspecialization,fill=HabitatType)) +
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=1.5) +
  stat_summary(fun=mean,geom="point",color="black",size=6,pch=21,stroke=2,aes(fill=HabitatType)) +
  ylab("Network Specialization (H2')") + xlab("Habitat Type") +
  ylim(0.1,0.5) +
  theme_classic() +
  scale_x_discrete(labels = c("Remnant\nForest",
                              "Restored\nForest")) +
  theme(axis.text.x=element_text(colour="black",size=14)) +
  theme(text=element_text(colour="black",size=14)) + 
  theme(axis.text.y=element_text(colour="black",size=14)) +
  theme(axis.title.y=element_text(colour="black",size=16, margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(colour="black", size=16, margin=margin(10,0,0,0))) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#e65262", "#39acd5")) 
  #stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  #ggtitle("A)") 

  ggsave("h2_plot.png", height = 6, width = 5)


#Observed welch t-test
ct_H2_welch <- t.test(subset(ct_H2_dat,HabitatType == "Restored")$networkspecialization,
                   subset(ct_H2_dat,HabitatType == "Remnant")$networkspecialization,
                   paired=FALSE, var.equal=FALSE)
ct_H2_welch


###########################
###   Linkage Density   ###
###########################

ct_link_dens <-lapply(haka_counts_bipartite, networklevel, index='linkage density')

# Mean number of interactions per species
ct_link_dens <- data.frame(unlist(ct_link_dens))
row.names(ct_link_dens) <- plot_names
ct_link_dens$plot <- plot_names
ct_link_dens$HabitatType <- rep(c("Restored","Remnant"),each=6)
names(ct_link_dens)[1] <- "LinkageDensity"



ct_link_dens_plot <- ggplot(ct_link_dens,aes(x=HabitatType,y=LinkageDensity)) +   
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=1.5) +
  stat_summary(fun=mean,geom="point",color="black",size=6,pch=21,stroke=2,aes(fill=HabitatType)) +
  ylab("Linkage Density") + xlab("Habitat Type") +
  scale_y_continuous(limits = c(7,15), breaks = seq(7, 15, by = 2)) +
  #coord_cartesian(ylim = c(7, 15)) +
  #ylim(8.5, 16) +
  theme_classic() +
  scale_x_discrete(labels = c("Remnant\nForest",
                              "Restored\nForest")) +
  theme(axis.text.x=element_text(colour="black",size=16,vjust=-1)) +
  theme(text=element_text(colour="black",size=16)) + 
  theme(axis.text.y=element_text(colour="black",size=16)) +
  theme(axis.title.y=element_text(colour="black",size=16)) +
  theme(axis.title.x=element_text(colour="white")) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#eed67c", "#8a9ee5")) +
  #stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  ggtitle("B)") 


#Observed welch t-test
ct_link_dens_welch <- t.test(subset(ct_link_dens,HabitatType == "Restored")$LinkageDensity,
                          subset(ct_link_dens,HabitatType == "Remnant")$LinkageDensity,
                          paired=FALSE, var.equal=FALSE)
ct_link_dens_welch


#######################
###   Connectance   ###
#######################
# ct_connectance <- lapply(haka_counts_bipartite, networklevel, index='connectance')
# 
# ct_con_dat <- data.frame(unlist(ct_connectance))
# row.names(ct_con_dat) <- plot_names
# ct_con_dat$plot <- plot_names
# ct_con_dat$HabitatType <- rep(c("Restored","Remnant"),each=6)
# names(ct_con_dat)[1] <- "connectance"
# 
# 
# ct_conn_plot <- ggplot(ct_con_dat,aes(x=HabitatType,y=connectance,fill=HabitatType)) +
#   stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",width=0.1,size=2) +
#   stat_summary(fun=mean,geom="point",size=6,shape=21,stroke=2) +
#   ylab("Network Connectance") + xlab("Habitat Type") +
#   ylim(0.2,0.5) +
#   xlab("") +
#   scale_x_discrete(labels = c("Remnant\nForest",
#                               "Restored\nForest")) +
#   theme(axis.text.x=element_text(colour="black",size=18,vjust=-1)) +
#   theme(text=element_text(colour="black",size=16)) + 
#   theme(axis.text.y=element_text(colour="black",size=16)) +
#   theme(axis.title.y=element_text(colour="black",size=20)) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#         panel.background=element_blank()) +
#   theme(legend.position="none") +
#   scale_fill_manual(values=c("#eed67c", "#8a9ee5")) +
#   #stat_compare_means(size=10,label="p.signif",label.x=1.5) +
#   ggtitle("C)") 
# 
# #Observed welch t-test
# ct_connectance_welch <- t.test(subset(ct_con_dat,HabitatType == "Restored")$connectance,
#                             subset(ct_con_dat,HabitatType == "Remnant")$connectance,
#                             paired=FALSE, var.equal=FALSE)
# ct_connectance_welch


### combine plots
all_plots <- ggarrange(ct_H2_plot, 
                       ct_link_dens_plot, 
                       #conn_plot, 
                       ncol = 2, 
                       nrow = 1) +
  theme(plot.margin = margin(0,1,0.5,0, "cm"))

plots_labeled <- annotate_figure(all_plots,
                                 bottom = text_grob("Habitat Type", 
                                                    size = 19))
ggsave("bipartite_plots.png", width = 10, height = 7)


### putting all metric results together

summarize_metric <- function (my_data, t_results) {
  
  remnant_mean <- mean(subset(my_data, HabitatType=="Remnant")[[1]])
  remnant_sd <- sd(subset(my_data, HabitatType=="Remnant")[[1]])
  
  restored_mean <- mean(subset(my_data, HabitatType=="Restored")[[1]])
  restored_sd <- sd(subset(my_data, HabitatType=="Restored")[[1]])
  
  ttest_pval <- t_results$p.value
  
  results <- c(remnant_mean, remnant_sd, restored_mean, restored_sd, ttest_pval)
  
  names(results) <- c("Remnant Mean",
                      "Remnant SD",
                      "Restored Mean",
                      "Restored SD",
                      "P")
  
  return(results)
  
}

ct_net_spec_results <- summarize_metric(my_data = ct_H2_dat,
                                     t_results = ct_H2_welch)

ct_link_dens_results <- summarize_metric(my_data = ct_link_dens,
                                      t_results = ct_link_dens_welch)

# conn_results <- summarize_metric(my_data = con_dat,
#                                  t_results = connectance_welch)

ct_bip_results <- t(data.frame(ct_net_spec_results, ct_link_dens_results))
ct_bip_results <- as.data.frame(ct_bip_results)


### P-value correction
library(multtest)

set.seed(20210917)

raw_p_bip = ct_bip_results$P

adj_p_bip = mt.rawp2adjp(raw_p_bip, na.rm = T)

unordered_p_bip = adj_p_bip$adjp[, c("rawp", "BH")]

ordered_p_bip = unordered_p_bip[,2][order(adj_p_bip$index)]

ct_bip_results$Adj_P <- ordered_p_bip

ct_bip_results$Index <- c("H2", "Linkage Density")

ct_bip_results <- select(ct_bip_results, 7, 1:6)

row.names(ct_bip_results) <- NULL

write.csv(ct_bip_results, "output/bipartite_metric_results_counts.csv")


#########################################
###   Host specialization on AMF (d') ###
#########################################

ct_dprime_vals <- list()

for (a_plot in plot_names) {
  
  host_spec <- dfun(t(haka_counts_bipartite[[a_plot]]))
  
  ct_dprime_vals[[a_plot]] <- host_spec$dprime
  
}



ct_host_spec_data <-as.data.frame(ct_dprime_vals)
ct_host_spec_df <- rownames_to_column(ct_host_spec_data, "Host")
ct_host_spec_df <- melt(ct_host_spec_df, id.vars = "Host")
ct_host_spec_df <- ct_host_spec_df %>% 
  dplyr::mutate(HabitatType = ifelse(grepl("AK", variable), "Restored", "Remnant"))
names(ct_host_spec_df)[2:3] <- c("plot", "dprime")

# quick test for if whole model is significant
# compare host d' values from restored to remnant
ct_dprime_test <- t.test(subset(ct_host_spec_df,HabitatType == "Restored")$dprime,
                         subset(ct_host_spec_df,HabitatType == "Remnant")$dprime,
                               paired=FALSE, var.equal=FALSE)
ct_dprime_test

# anova version
a <- aov(formula = dprime ~ HabitatType, data = ct_host_spec_df)
summary(a)

# comparing d' across hosts in general
h <- aov(formula = dprime ~ Host, data = ct_host_spec_df)
summary(h)
b = TukeyHSD(h)

# within habitats
ro_dprime <- ct_host_spec_df[ct_host_spec_df$HabitatType=="Remnant",]
ak_dprime <- ct_host_spec_df[ct_host_spec_df$HabitatType=="Restored",]

ro_d_aov <- aov(formula = dprime ~ Host, data = ro_dprime)
ak_d_aov <- aov(formula = dprime ~ Host, data = ak_dprime)

ro_d_tuk <- TukeyHSD(ro_d_aov)
ak_d_tuk <- TukeyHSD(ak_d_aov)



# by host and habitat
tuk_dprime_df <- ct_host_spec_df

tuk_dprime_df$hab_host <- paste(substr(tuk_dprime_df$plot, 1, 2),
                                tuk_dprime_df$Host,
                                sep = "_")
hab_host_aov <- aov(formula = dprime ~ hab_host, data = tuk_dprime_df)
hab_host_tuk <- TukeyHSD(hab_host_aov)

tuk_df <- as.data.frame(hab_host_tuk$hab_host)
tuk_df <- rownames_to_column(tuk_df, "Habitat_Host")
tuk_df <- tuk_df[,c(1,2,5)]
write.csv(tuk_df, "output/tukey_outputs.csv", row.names = F)


ct_host_spec_t <- list()

#full_host_names <- unique(host_spec_df$Host)

host_names_no_soil <- host_names[1:6]

for (a_host in host_names_no_soil) {
  
  remnant_vals <- ct_host_spec_df$dprime[ct_host_spec_df$Host==a_host & 
                                          ct_host_spec_df$HabitatType=="Remnant"]
  
  restored_vals <- ct_host_spec_df$dprime[ct_host_spec_df$Host==a_host & 
                                            ct_host_spec_df$HabitatType=="Restored"]
  
  t <- t.test(remnant_vals, restored_vals)
  
  ct_host_spec_t[[a_host]] <- t
  
}


ct_d_df <- as.data.frame(host_names_no_soil)
names(ct_d_df) <- "Host"

ct_d_list <- list()

for(a_host in host_names_no_soil) {
  
  remnant_mean <- mean(subset(ct_host_spec_df, Host == a_host & HabitatType == "Remnant")$dprime)
  remnant_sd <- sd(subset(ct_host_spec_df, Host == a_host & HabitatType == "Remnant")$dprime)
  
  restored_mean <- mean(subset(ct_host_spec_df, Host == a_host & HabitatType == "Restored")$dprime)
  restored_sd <- sd(subset(ct_host_spec_df, Host == a_host & HabitatType == "Restored")$dprime)
  
  d_pval <- ct_host_spec_t[[a_host]]$p.value
  
  d_vals <- c(remnant_mean, remnant_sd, restored_mean, restored_sd, d_pval)
  
  names(d_vals) <- c("Remnant Mean",
                     "Remnant SD",
                     "Restored Mean",
                     "Restored SD",
                     "P")
  
  ct_d_list[[a_host]] <- d_vals
  
}

ct_dprime_table <- do.call(rbind, ct_d_list)
ct_dprime_table <- as.data.frame(ct_dprime_table)
ct_dprime_table <- rownames_to_column(ct_dprime_table, "Host")

### P-value correction

set.seed(20210917)

ct_raw_p_values = ct_dprime_table$P

ct_adj_p_values = mt.rawp2adjp(ct_raw_p_values, na.rm = T)

ct_unordered_p = ct_adj_p_values$adjp[, c("rawp", "BH")]

ct_ordered_p = ct_unordered_p[,2][order(ct_adj_p_values$index)]

ct_dprime_table$Adj_P <- ct_ordered_p

write.csv(ct_dprime_table, "output/dprime_table.csv")



ct_host_dfun_p <- ggplot(ct_host_spec_df,aes(x=HabitatType,y=dprime,fill=factor(Host, levels=host_names))) + 
  geom_boxplot(size=0.75, outlier.alpha = 0) +
  geom_point(pch=21,position=position_jitterdodge(jitter.width=0.1),size=3) +
  labs(x = "Habitat Type", y = "Host Specialization on AM fungi (d')", fill="Host") +
  ylim(0,0.8) +
  #scale_fill_manual(values=c("#D73027","#FC8D59","#FEE090","#008000","#E0F3F8","#008B8B","#4575B4")) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c("Remnant\nForest",
                              "Restored\nForest")) +
  theme(text=element_text(colour="black",size=14)) + 
  theme(legend.key = element_rect(fill="white")) +
  theme(axis.text.x=element_text(colour="black",size=12)) +
  theme(axis.text.y=element_text(colour="black",size=12)) +
  theme(axis.title.x=element_text(colour="black",size=14, margin = margin(10,0,0,0))) +
  theme(axis.title.y=element_text(colour="black",size=14, margin = margin(0,10,0,0))) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.text = element_text(face="italic",size=12)) 

ggsave("dprime_plot_using_counts.png", width = 7.5, height = 6)



