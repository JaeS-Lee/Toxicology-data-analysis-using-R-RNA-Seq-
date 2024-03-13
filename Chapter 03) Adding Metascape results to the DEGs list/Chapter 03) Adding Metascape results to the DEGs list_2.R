library(tidyverse)
library(readxl)
library(writexl)

## Load DEGs excel files that we made in Chapter 1

A_1_control_DEGs <- readxl::read_excel('./DEGs/01_2.A_1_control_DEGs.xlsx')
A_1_control_DEGs <- A_1_control_DEGs %>% dplyr::rename(MyList = Ensembl_ID)
head(A_1_control_DEGs)

A_10_control_DEGs <- readxl::read_excel('./DEGs/02_2.A_10_control_DEGs.xlsx')
A_10_control_DEGs <- A_10_control_DEGs %>% dplyr::rename(MyList = Ensembl_ID)
head(A_10_control_DEGs)


## Load Metascape excel files that we downloaded in Chapter 2

A_1_control_metascape <- readxl::read_excel('./Metascape results/01_2.A_1_control_DEGs/metascape_result_A_1_control.xlsx', sheet = 1)
A_1_control_metascape <- A_1_control_metascape %>% dplyr::select(MyList, `Gene Symbol`, Description) %>% dplyr::filter(`Gene Symbol` != 'None')
A_1_control_metascape <- A_1_control_metascape %>% dplyr::rename(Gene_Symbol = `Gene Symbol`)
head(A_1_control_metascape)

A_10_control_metascape <- readxl::read_excel('./Metascape results/02_2.A_10_control_DEGs/metascape_result_A_10_control.xlsx', sheet = 1)
A_10_control_metascape <- A_10_control_metascape %>% dplyr::select(MyList, `Gene Symbol`, Description) %>% dplyr::filter(`Gene Symbol` != 'None')
A_10_control_metascape <- A_10_control_metascape %>% dplyr::rename(Gene_Symbol = `Gene Symbol`)
head(A_10_control_metascape)


## Merging two separated data

A_1_control_merged <- merge(A_1_control_DEGs, A_1_control_metascape, by = c('MyList'))
A_1_control_merged <- A_1_control_merged %>% select(MyList, Gene_Symbol, logFC, FDR) %>% arrange(desc(logFC))
#writexl::write_xlsx(A_1_control_merged, "./Output/01-2.A_1_control_merged.xlsx")

A_10_control_merged <- merge(A_10_control_DEGs, A_10_control_metascape, by = c('MyList'))
A_10_control_merged <- A_10_control_merged %>% select(MyList, Gene_Symbol, logFC, FDR) %>% arrange(desc(logFC))
#writexl::write_xlsx(A_10_control_merged, "./Output/02-2.A_10_control_merged.xlsx")


## Merging all data in a chart (Outer join)

A_1_control_merged_for_outer <- A_1_control_merged %>% 
  dplyr::rename(A_1_logFC = logFC, A_1_FDR = FDR)

A_10_control_merged_for_outer <- A_10_control_merged %>% 
  dplyr::rename(A_10_logFC = logFC, A_10_FDR = FDR)


outer1 <- merge(A_1_control_merged_for_outer, A_10_control_merged_for_outer, by=c('MyList', 'Gene_Symbol'), all = TRUE)
outer1 <- outer1 %>% dplyr::arrange(desc(A_1_logFC), desc(A_10_logFC))


# LogFC saver

outer1_logFC <- outer1 %>% dplyr::select(MyList, Gene_Symbol, A_1_logFC, A_10_logFC)
View(outer1_logFC)
#writexl::write_xlsx(outer1_logFC, "./Output/03.All_outer_merged_logFC.xlsx")