library(tidyverse)
library(data.table)
library(writexl)


## Load Metascape excel files

A_1_control_metascape <- readxl::read_excel('./Metascape results/01_2.A_1_control_DEGs/metascape_result_A_1_control.xlsx', sheet = 2)
A_10_control_metascape <- readxl::read_excel('./Metascape results/02_2.A_10_control_DEGs/metascape_result_A_10_control.xlsx', sheet = 2)

result_list <- list()

result_list[['A_1_control']] <- A_1_control_metascape
result_list[['A_10_control']] <- A_10_control_metascape

## Select GroupID = Summary, add group names, description, and Genes

selected_list = list()

for(j in names(result_list)){
  temp2 <- c()
  # Remove Summary parts because the information about members is necessary
  temp2 <- dplyr::filter(result_list[[j]], str_detect(GroupID,'Summary') != T) 
  temp2 <- temp2 %>% dplyr::mutate(Group = j)
  temp2 <- temp2 %>% dplyr::select(Group, Category, Description, Genes)
  
  selected_list[[paste(j, c("summary_deleted"), sep = '_')]] <- temp2
}

selected_list


## Show all the data in 'selected_list' with head function

##########################################################################################################################
## Summary of enrichment result by chemicals
##########################################################################################################################

# A_1_control

names(selected_list)

total_relocation_A_1_control <- data.table() # blank data.table
preprocessed_total_A_1_control <- selected_list$A_1_control_summary_deleted
head(as_tibble(preprocessed_total_A_1_control))

total_relocation_row_num_A_1_control <- nrow(preprocessed_total_A_1_control) # total number of rows

for (i in 1:total_relocation_row_num_A_1_control){
  
  Group_index_A_1_control <- as.character(preprocessed_total_A_1_control[i, 1])
  Category_index_A_1_control <- as.character(preprocessed_total_A_1_control[i, 2])
  Description_index_A_1_control <- as.character(preprocessed_total_A_1_control[i, 3])
  Genes_index_A_1_control <- as.character(preprocessed_total_A_1_control[i, 4])
  Genes_index_split_temp_A_1_control <- data.frame(strsplit(Genes_index_A_1_control, split = ','))
  relocation_temp_A_1_control <- data.frame(cbind(Group_index_A_1_control, Category_index_A_1_control, Description_index_A_1_control, 
                                                  Genes_index_split_temp_A_1_control))
  names(relocation_temp_A_1_control) <- c("Group", 'Category', "Description", "Genes")
  total_relocation_A_1_control  <- rbind(total_relocation_A_1_control, relocation_temp_A_1_control)
}

total_relocation_A_1_control # 

total_result_A_1_control <- total_relocation_A_1_control
total_result_A_1_control <- total_result_A_1_control %>% dplyr::arrange(Category)
total_result_A_1_control

#writexl::write_xlsx(total_result_A_1_control,"./Output/01-1.A_1_control_Summary_Enrichment_genes.xlsx")

# Reference : https://rfriend.tistory.com/238

##########################################################################################################################

# A_10_control

names(selected_list)

total_relocation_A_10_control <- data.table() # blank data.table
preprocessed_total_A_10_control <- selected_list$A_10_control_summary_deleted
head(as_tibble(preprocessed_total_A_10_control))

total_relocation_row_num_A_10_control <- nrow(preprocessed_total_A_10_control) # total number of rows

for (i in 1:total_relocation_row_num_A_10_control){
  
  Group_index_A_10_control <- as.character(preprocessed_total_A_10_control[i, 1])
  Category_index_A_10_control <- as.character(preprocessed_total_A_10_control[i, 2])
  Description_index_A_10_control <- as.character(preprocessed_total_A_10_control[i, 3])
  Genes_index_A_10_control <- as.character(preprocessed_total_A_10_control[i, 4])
  Genes_index_split_temp_A_10_control <- data.frame(strsplit(Genes_index_A_10_control, split = ','))
  relocation_temp_A_10_control <- data.frame(cbind(Group_index_A_10_control, Category_index_A_10_control, Description_index_A_10_control, 
                                                  Genes_index_split_temp_A_10_control))
  names(relocation_temp_A_10_control) <- c("Group", 'Category', "Description", "Genes")
  total_relocation_A_10_control  <- rbind(total_relocation_A_10_control, relocation_temp_A_10_control)
}

total_relocation_A_10_control # 

total_result_A_10_control <- total_relocation_A_10_control
total_result_A_10_control <- total_result_A_10_control %>% dplyr::arrange(Category)
total_result_A_10_control

#writexl::write_xlsx(total_result_A_10_control,"./Output/02-1.A_10_control_Summary_Enrichment_genes.xlsx")

##########################################################################################################################