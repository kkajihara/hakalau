
otu_table_prep <- function (a_host, physeq_list, ro_grouping, ak_grouping) {
  
  host_otu_table <- otu_table(physeq_list[[a_host]])
  
  host_otu_table <- as.data.frame(host_otu_table)
  
  host_otu_table$Host <- rep(a_host, nrow(host_otu_table))
  
  host_otu_table <- rownames_to_column(host_otu_table, "ASV")
  
  ak_sum_for_rel_abun <- sum(host_otu_table[grepl("AK", colnames(host_otu_table))])
  
  host_otu_table[,(grepl("AK", names(host_otu_table)))] <-
    apply(host_otu_table[,(grepl("AK", names(host_otu_table)))], 
          2,
          function(x) (x/ak_sum_for_rel_abun)*100)
  
  ro_sum_for_rel_abun <- sum(host_otu_table[grepl("RO", colnames(host_otu_table))])
  
  host_otu_table[,(grepl("RO", names(host_otu_table)))] <-
    apply(host_otu_table[,(grepl("RO", names(host_otu_table)))],
          2,
          function(x) (x/ro_sum_for_rel_abun)*100)
  
  host_otu_table$Family <- full_tax_table_df$Family[match(
    host_otu_table$ASV, full_tax_table_df$ASV)]
  
  ak_host_table <- host_otu_table[,!(grepl("RO", names(host_otu_table)))]
  
  ak_unique_host_table <- ak_host_table[ak_host_table$ASV %in% ak_grouping[[a_host]]$ASV,]
  
  ak_unique_host_table <- 
    ak_unique_host_table[, colSums(ak_unique_host_table != 0) > 0]
  
  
  ak_unique_host_table <- melt(ak_unique_host_table)
  
  colnames(ak_unique_host_table)[c(4,5)] <- c("SampleID", "rel_abun")
  
  ak_unique_host_table <- ak_unique_host_table[ak_unique_host_table$rel_abun > 0,]
  
  
  ro_host_table <- host_otu_table[,!(grepl("AK", names(host_otu_table)))]
  
  ro_unique_host_table <- ro_host_table[ro_host_table$ASV %in% ro_grouping[[a_host]]$ASV,]
  
  ro_unique_host_table <- melt(ro_unique_host_table)
  
  colnames(ro_unique_host_table)[c(4,5)] <- c("SampleID", "rel_abun")
  
  ro_unique_host_table <- ro_unique_host_table[ro_unique_host_table$rel_abun > 0,]
  
  full_host_table <- rbind(ak_unique_host_table, ro_unique_host_table)
  
  full_host_table$Habitat <- substr(full_host_table$SampleID, 
                                    start = 1,
                                    stop = 2)
  
  
  return(full_host_table)
  
  
}




anova_prep <- function (a_host, relabun_df) {
  
  host_df <- relabun_df[[a_host]] %>%
    group_by(SampleID, Habitat, Family, Host) %>%
    summarise(rel_abun=sum(rel_abun))
  
  long_host_df <- host_df %>%
    dcast(SampleID + Habitat ~ Family, value.var = "rel_abun")
  
  return(long_host_df)
  
}


# function to run linear models on families' relative abundance within host data tables
# Habitat is the independent variable
run_lm <- function(family_name, lm_abun) {
  
  while (family_name %in% names(lm_abun)) {
    
    formula_call <- as.formula(paste0(family_name, " ~ Habitat"))
    
    lm_out <- lm(formula_call,
                 data = lm_abun)
    return(lm_out)
  }
  
}



summarize_lm <- function(family, lm_out) {
  
  while (family %in% names(lm_out)) {
    
    #for (a_host in host_names) {
    # summarize model output
    lm_sum <- summary(aov(lm_out[[family]]))
    
    #}
    
    summary_out <- 
      data.frame(
        p_val  = lm_sum[[1]][["Pr(>F)"]][[1]]
      )
    
    return(summary_out)
    
  }
  
}




null_remove <- function (a_host, a_lm) {
  
  no_null <- a_lm[[a_host]][!sapply(a_lm[[a_host]], is.null)]
  
  return(no_null)
  
}


null_to_na <- function (a_host, a_aov) {
  
  nullstoNA <- a_aov[[a_host]]
  
  nullstoNA[sapply(nullstoNA, is.null)] <- NA
  
  return(nullstoNA)
  
}


anova_tbl_format <- function (pval_table) {
  
  # round each value to 4 digits
  p_table <- apply(pval_table, c(1,2), function(x) round(x, digits = 4))
  
  # transpose so hosts are columns, families are rows
  p_table <- t(p_table)
  
  return(p_table)
  
}







