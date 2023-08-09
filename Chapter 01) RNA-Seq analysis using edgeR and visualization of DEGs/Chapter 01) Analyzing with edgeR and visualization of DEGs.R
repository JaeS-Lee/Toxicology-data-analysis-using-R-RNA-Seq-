library(tidyverse)
library(edgeR)
library(EnhancedVolcano)
library(ggvenn)
library(writexl)

CountData <- read.csv("./dummy_dataset.csv", row.names = 1)
head(CountData)

group <- factor(c('Control', 'Control', 'A_1uM', 'A_1uM', 'A_10uM', 'A_10uM'),
                levels = unique(c('Control', 'Control', 'A_1uM', 'A_1uM', 'A_10uM', 'A_10uM')))

# Reference : https://stackoverflow.com/questions/53195124/error-factor-level-2-is-duplicated-when-not-really-duplicated

design <- model.matrix(~ group)
design

d <- DGEList(counts = CountData, group = group)

keep_data <- filterByExpr(d) # Filtering right after DGEList
table(keep_data)

d <- d[keep_data, ,keep.lib.sizes = FALSE]
d

# Reference : https://heavywatal.github.io/rstats/edger.html

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

fit <- glmFit(d, design)

###################################################################################################################

## A_1 vs DMSO

lrt_A_1 <- glmLRT(fit, coef = c(2))
topTags(lrt_A_1)

total_table_A_1_control <- as.data.frame(topTags(lrt_A_1, n = nrow(lrt_A_1)))
total_table_A_1_control <- total_table_A_1_control %>% tibble::rownames_to_column(var = 'Ensembl_ID')
#writexl::write_xlsx(total_table_A_1_control, "./Results/01.Raw/01_1.A_1_control_Raw.xlsx")

total_table_A_1_control_DEGs <- total_table_A_1_control %>% dplyr::filter(abs(logFC) > 1 & FDR < 0.05)
head(total_table_A_1_control_DEGs);nrow(total_table_A_1_control_DEGs)
#writexl::write_xlsx(total_table_A_1_control_DEGs, "./Results/02.DEGs/01_2.A_1_control_DEGs.xlsx")


# Volcano plot A_1

max(total_table_A_1_control$logFC);min(total_table_A_1_control$logFC)
-log10(max(total_table_A_1_control$FDR));-log10(min(total_table_A_1_control$FDR))

keyvals_A_1 <- ifelse(
  total_table_A_1_control$logFC < -1 & total_table_A_1_control$FDR < 0.05, 'red',
  ifelse(total_table_A_1_control$logFC > 1 & total_table_A_1_control$FDR < 0.05, 'forestgreen', 'black'))

keyvals_A_1[is.na(keyvals_A_1)] <- 'black'
names(keyvals_A_1)[keyvals_A_1 == 'forestgreen'] <- 'Upregulated'
names(keyvals_A_1)[keyvals_A_1 == 'black'] <- 'Non-DEGs'
names(keyvals_A_1)[keyvals_A_1 == 'red'] <- 'Downregulated'

EnhancedVolcano(total_table_A_1_control,
                lab = NA,
                x = 'logFC',
                y = 'FDR',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote('-'~Log[10]~ 'FDR'),
                title = 'A_1 volcano plot',
                ylim = c(0, 299),
                xlim = c(-10, 10),
                colCustom = keyvals_A_1)

###################################################################################################################

## A_10 vs DMSO

lrt_A_10 <- glmLRT(fit, coef = c(3))
topTags(lrt_A_10)

total_table_A_10_control <- as.data.frame(topTags(lrt_A_10, n = nrow(lrt_A_10)))
total_table_A_10_control <- total_table_A_10_control %>% tibble::rownames_to_column(var = 'Ensembl_ID')
#writexl::write_xlsx(total_table_A_10_control, "./Results/01.Raw/02_1.A_10_control_Raw.xlsx")

total_table_A_10_control_DEGs <- total_table_A_10_control %>% dplyr::filter(abs(logFC) > 1 & FDR < 0.05)
head(total_table_A_10_control_DEGs);nrow(total_table_A_10_control_DEGs)
#writexl::write_xlsx(total_table_A_10_control_DEGs, "./Results/01.Raw/02_2.A_10_control_DEGs.xlsx")


# Volcano plot A_10

max(total_table_A_10_control$logFC);min(total_table_A_10_control$logFC)
-log10(max(total_table_A_10_control$FDR));-log10(min(total_table_A_10_control$FDR))

keyvals_A_10 <- ifelse(
  total_table_A_10_control$logFC < -1 & total_table_A_10_control$FDR < 0.05, 'red',
  ifelse(total_table_A_10_control$logFC > 1 & total_table_A_10_control$FDR < 0.05, 'forestgreen', 'black'))

keyvals_A_10[is.na(keyvals_A_10)] <- 'black'
names(keyvals_A_10)[keyvals_A_10 == 'forestgreen'] <- 'Upregulated'
names(keyvals_A_10)[keyvals_A_10 == 'black'] <- 'Non-DEGs'
names(keyvals_A_10)[keyvals_A_10 == 'red'] <- 'Downregulated'

EnhancedVolcano(total_table_A_10_control,
                lab = NA,
                x = 'logFC',
                y = 'FDR',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote('-'~Log[10]~ 'FDR'),
                title = 'A_10 volcano plot',
                ylim = c(0, 299),
                xlim = c(-10, 10),
                colCustom = keyvals_A_10)

# One or more p-values is 0. Converting to 10^-1 * current lowest non-zero p-value... 

###################################################################################################################

## Venn diagram (Up and Down)

# Venndiagram from DEGs 

result_list_Venn <- list()

A_1_DEGs_Venn <- total_table_A_1_control_DEGs
head(A_1_DEGs_Venn)

A_10_DEGs_Venn <- total_table_A_10_control_DEGs
head(A_10_DEGs_Venn)

result_list_Venn[['A_1_DEGs_Venn']] <- A_1_DEGs_Venn
result_list_Venn[['A_10_DEGs_Venn']] <- A_10_DEGs_Venn

head(result_list_Venn)

#############################################################################################

## A_1 DEGs for Venndiagram

A_1_DEGs_Venn_df <- result_list_Venn[[1]]
sum(abs(A_1_DEGs_Venn_df$logFC) > 1) == nrow(A_1_DEGs_Venn_df) # True
sum(abs(A_1_DEGs_Venn_df$FDR) < 0.05) == nrow(A_1_DEGs_Venn_df) # True
A_1_DEGs_Venn_df <- A_1_DEGs_Venn_df %>% dplyr::select(Ensembl_ID, logFC)

# A_1 upregulated DEGs list
A_1_DEGs_Venn_Up_vector <- A_1_DEGs_Venn_df %>% dplyr::filter(logFC > 1)
A_1_DEGs_Venn_Up_vector <- A_1_DEGs_Venn_Up_vector %>% dplyr::select(Ensembl_ID)
A_1_DEGs_Venn_Up_vector <- A_1_DEGs_Venn_Up_vector[[1]]

# A_1 Downregulated DEGs list
A_1_DEGs_Venn_Down_vector <- A_1_DEGs_Venn_df %>% dplyr::filter(logFC < -1)
A_1_DEGs_Venn_Down_vector <- A_1_DEGs_Venn_Down_vector %>% dplyr::select(Ensembl_ID)
A_1_DEGs_Venn_Down_vector <- A_1_DEGs_Venn_Down_vector[[1]]

########################################

## A_10 DEGs for Venndiagram

A_10_DEGs_Venn_df <- result_list_Venn[[2]]
sum(abs(A_10_DEGs_Venn_df$logFC) > 1) == nrow(A_10_DEGs_Venn_df) # True
sum(abs(A_10_DEGs_Venn_df$FDR) < 0.05) == nrow(A_10_DEGs_Venn_df) # True
A_10_DEGs_Venn_df <- A_10_DEGs_Venn_df %>% dplyr::select(Ensembl_ID, logFC)

# A_10 upregulated DEGs list

A_10_DEGs_Venn_Up_vector <- A_10_DEGs_Venn_df %>% dplyr::filter(logFC > 1)
A_10_DEGs_Venn_Up_vector <- A_10_DEGs_Venn_Up_vector %>% dplyr::select(Ensembl_ID)
A_10_DEGs_Venn_Up_vector <- A_10_DEGs_Venn_Up_vector[[1]]

# A_10 Downregulated DEGs list

A_10_DEGs_Venn_Down_vector <- A_10_DEGs_Venn_df %>% dplyr::filter(logFC < -1)
A_10_DEGs_Venn_Down_vector <- A_10_DEGs_Venn_Down_vector %>% dplyr::select(Ensembl_ID)
A_10_DEGs_Venn_Down_vector <- A_10_DEGs_Venn_Down_vector[[1]]

#############################################################################################

## Venndiagram (Figures) 500 x 500

# Upregulation (A_1 & A_10)

All_up <- list('A_1' = A_1_DEGs_Venn_Up_vector,
               'A_10' = A_10_DEGs_Venn_Up_vector)

ggvenn(All_up, 
       show_percentage=FALSE, 
       fill_color = c('green','forestgreen'),
       fill_alpha = 0.7,
       stroke_color = NA,
       text_size = 10)


# Downregulation (A_1 & A_10)

All_down <- list('A_1' = A_1_DEGs_Venn_Down_vector,
                 'A_10' = A_10_DEGs_Venn_Down_vector)

ggvenn(All_down, 
       show_percentage=FALSE, 
       fill_color = c('orange','red'),
       fill_alpha = 0.7,
       stroke_color = NA,
       text_size = 10)
