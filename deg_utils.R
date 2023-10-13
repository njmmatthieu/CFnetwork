deg_in_common <- function(deg_list, nb_in_common=5, studies_order=NA) {
  
  all_UP_deg <- lapply(deg_list, function(deg_df){
    return(deg_df[which(deg_df$diffexpressed=="UP"),"Symbol"])
  })
  all_DOWN_deg <- lapply(deg_list, function(deg_df){
    return(deg_df[which(deg_df$diffexpressed=="DOWN"),"Symbol"])
  })
  
  diff_deg <- c(unlist(all_UP_deg),
                         unlist(all_DOWN_deg))
  
  deg_occurrence <- data.frame(table(diff_deg))
  colnames(deg_occurrence) <- c("Symbol", "Occurrence")
  
  intersect_deg <- as.character(deg_occurrence[which(deg_occurrence$Occurrence>=nb_in_common),"Symbol"])
  
  diff_intersections_summary <- data.frame(t(sapply(intersect_deg, function(symbol){
    up_deg <- unlist(lapply(all_UP_deg, function(up_deg) as.integer(symbol %in% up_deg)))
    down_deg <- unlist(lapply(all_DOWN_deg, function(down_deg) as.integer(symbol %in% down_deg)))
    
    summary_deg <- up_deg - down_deg
    return(summary_deg)
  })))

  colnames(diff_intersections_summary) <- names(deg_list)
  rownames(diff_intersections_summary) <- intersect_deg
  
  if (!is.na(studies_order)){
    diff_intersections_summary <- diff_intersections_summary[,studies_order]
  }
  
  return(diff_intersections_summary)
  
}

deg_pairwise_in_common <- function(diff_intersections_summary) {
  
    deg_pairwise_in_common <- data.frame()
    for (i_study_A in 1:(ncol(diff_intersections_summary)-1)){
      study_A <- colnames(diff_intersections_summary)[i_study_A]
      for (i_study_B in (i_study_A+1):ncol(diff_intersections_summary)){
        study_B <- colnames(diff_intersections_summary[i_study_B])
        
        diff_pairwise <- diff_intersections_summary[,c(study_A,study_B)]
        diff_pairwise <- diff_pairwise[which(apply(diff_pairwise, 1, function(t) !all(t==0))),]
        
        common_UP <- length(which(rowSums(diff_pairwise[,c(study_A,study_B)])==2))
        common_BOTH <- length(which(rowSums(diff_pairwise[,c(study_A,study_B)])==0))
        common_DOWN <- length(which(rowSums(diff_pairwise[,c(study_A,study_B)])==-2))
        deg_pairwise <- data.frame(study_A = study_A,
                                   study_B = study_B,
                                   pair = paste(study_A, study_B, sep = "_"),
                                   UP = common_UP,
                                   BOTH = common_BOTH,
                                   DOWN = common_DOWN,
                                   ALL = common_UP + common_BOTH + common_DOWN)
        deg_pairwise_in_common <- rbind(deg_pairwise_in_common,
                                        deg_pairwise)
      }
    }
    
    deg_pairwise_in_common_long <- melt(deg_pairwise_in_common,id.vars = c("study_A","study_B","ALL", "pair"))
    deg_pairwise_in_common_long$prop <- deg_pairwise_in_common_long$value / deg_pairwise_in_common_long$ALL
    deg_pairwise_in_common_long$prop[which(is.na(deg_pairwise_in_common_long$prop))]<-0
    
    
    deg_pairwise_in_common_long <- deg_pairwise_in_common_long %>% group_by(pair) %>%
      mutate(cp1=c(0,head(cumsum(prop),-1)),
             cp2=cumsum(prop))
    
    deg_pairwise_in_common_long$study_A <- factor(deg_pairwise_in_common_long$study_A,
                                                  levels = studies_names)
    deg_pairwise_in_common_long$study_B <- factor(deg_pairwise_in_common_long$study_B,
                                                  levels = studies_names)
    
    deg_pairwise_in_common_long$log_all <- log(deg_pairwise_in_common_long$ALL)
    
    return(deg_pairwise_in_common_long)
  
}

