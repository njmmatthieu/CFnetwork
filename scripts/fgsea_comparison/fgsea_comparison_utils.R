# Keep from the matrix of NES of all pathways, the pathways DEP in at least 3 studies
fgsea_nes_diff_pathways_in_common <- function(fgsea_lists, 
                                              fgsea_combined_df, 
                                              nb_in_common=5, 
                                              studies_order=NA) {
  
  all_UP_pathways <- lapply(fgsea_lists, function(module_df){
    return(module_df[which(module_df$diff_activated=="UP"),"pathway"])
  })
  all_DOWN_pathways <- lapply(fgsea_lists, function(module_df){
    return(module_df[which(module_df$diff_activated=="DOWN"),"pathway"])
  })
  
  # print(all_UP_pathways)
  
  all_diff_pathways <- c(unlist(all_UP_pathways),
                         unlist(all_DOWN_pathways))
  
  
  
  diff_pathways_occurrence <- data.frame(table(all_diff_pathways))
  colnames(diff_pathways_occurrence) <- c("Pathway", "Occurrence")
  
  intersect_pathways <- as.character(diff_pathways_occurrence[which(diff_pathways_occurrence$Occurrence>=nb_in_common),"Pathway"])
  
  # print(length(intersect_pathways))
  
  fgsea_diff_pathways_combined_df <- fgsea_combined_df[intersect_pathways,]
  
  if (!is.na(studies_order)){
    fgsea_diff_pathways_combined_df <- fgsea_diff_pathways_combined_df[,studies_order]
  }
  
  return(fgsea_diff_pathways_combined_df)
  
}

# Keep from the matrix of NES of the pathways DEP in at least 3 studies:
# 1: if the DEP is significantly UP
# 0: not DEP in the particular study
# -1: if the DEP is significantly DOWN
fgsea_diff_pathways_in_common <- function(fgsea_lists, nb_in_common=5, studies_order=NA) {
  
  all_UP_pathways <- lapply(fgsea_lists, function(module_df){
    # print(module_df[which(module_df$diff_activated=="UP"),"pathway"])
    return(module_df[which(module_df$diff_activated=="UP"),"pathway"])
  })
  all_DOWN_pathways <- lapply(fgsea_lists, function(module_df){
    return(module_df[which(module_df$diff_activated=="DOWN"),"pathway"])
  })
  
  
  all_diff_pathways <- c(unlist(all_UP_pathways),
                         unlist(all_DOWN_pathways))
  
  diff_pathways_occurrence <- data.frame(table(all_diff_pathways))
  colnames(diff_pathways_occurrence) <- c("Pathway", "Occurrence")
  
  intersect_pathways <- as.character(diff_pathways_occurrence[which(diff_pathways_occurrence$Occurrence>=nb_in_common),"Pathway"])
  
  
  diff_intersections_summary <- sapply(intersect_pathways, function(pathway){
    up_pathways <- unlist(lapply(all_UP_pathways, function(pathways_up){
      return(as.integer(pathway %in% unlist(pathways_up)))
      }))
    down_pathways <- unlist(lapply(all_DOWN_pathways, function(pathways_down) {
      return(as.integer(pathway %in% unlist(pathways_down)))
    }))
    
    summary_pathways <- up_pathways - down_pathways
    return(summary_pathways)
  })
  
  diff_intersections_summary_df <- data.frame(t(diff_intersections_summary))
  colnames(diff_intersections_summary_df) <- names(fgsea_lists)
  rownames(diff_intersections_summary_df) <- intersect_pathways
  
  if (!is.na(studies_order)){
    diff_intersections_summary_df <- diff_intersections_summary_df[,studies_order]
  }
  
  return(diff_intersections_summary_df)
  
}

fgsea_dep_pairwise_in_common <- function(diff_intersections_summary) {
  
  dep_pairwise_in_common <- data.frame()
  for (i_study_A in 1:(ncol(diff_intersections_summary)-1)){
    study_A <- colnames(diff_intersections_summary)[i_study_A]
    for (i_study_B in (i_study_A+1):ncol(diff_intersections_summary)){
      study_B <- colnames(diff_intersections_summary[i_study_B])
      
      diff_pairwise <- diff_intersections_summary[,c(study_A,study_B)]
      diff_pairwise <- diff_pairwise[which(apply(diff_pairwise, 1, function(t) !all(t==0))),]
      
      common_UP <- length(which(rowSums(diff_pairwise[,c(study_A,study_B)])==2))
      common_BOTH <- length(which(rowSums(diff_pairwise[,c(study_A,study_B)])==0))
      common_DOWN <- length(which(rowSums(diff_pairwise[,c(study_A,study_B)])==-2))
      dep_pairwise <- data.frame(study_A = study_A,
                                 study_B = study_B,
                                 pair = paste(study_A, study_B, sep = "_"),
                                 UP = common_UP,
                                 BOTH = common_BOTH,
                                 DOWN = common_DOWN,
                                 ALL = common_UP + common_BOTH + common_DOWN)
      dep_pairwise_in_common <- rbind(dep_pairwise_in_common,
                                      dep_pairwise)
    }
  }
  
  dep_pairwise_in_common_long <- melt(dep_pairwise_in_common,id.vars = c("study_A","study_B","ALL", "pair"))
  dep_pairwise_in_common_long$prop <- dep_pairwise_in_common_long$value / dep_pairwise_in_common_long$ALL
  dep_pairwise_in_common_long$prop[which(is.na(dep_pairwise_in_common_long$prop))]<-0
  
  
  dep_pairwise_in_common_long <- dep_pairwise_in_common_long %>% group_by(pair) %>%
    mutate(cp1=c(0,head(cumsum(prop),-1)),
           cp2=cumsum(prop))
  
  dep_pairwise_in_common_long$study_A <- factor(dep_pairwise_in_common_long$study_A,
                                                levels = colnames(diff_intersections_summary))
  dep_pairwise_in_common_long$study_B <- factor(dep_pairwise_in_common_long$study_B,
                                                levels = colnames(diff_intersections_summary))
  
  dep_pairwise_in_common_long$log_all <- log(dep_pairwise_in_common_long$ALL)
  
  return(dep_pairwise_in_common_long)
  
}