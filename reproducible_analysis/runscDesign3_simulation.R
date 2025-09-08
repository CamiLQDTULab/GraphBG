rm(list = ls())
gc()  # Optional: Ask R to run garbage collection and free up memory

library(tidyverse)
library(viridis)
#install.packages("ggseg")
library(ggseg)
library(sf)
library(sp)
# n_genes = 1000, ADT = 5
# n_genes = 10000, ADT = 10
# n_cells = 2k, 5k
sampling_rate <- 7.9 # 20 = 2k cells, 7.9 (old: 8) for 5k cells, 3.95 (old: 3.9) for 10k cells, 40 for 1k cells
n_genes <- 100
data_name <- "CBMCs"


p <- ggplot() +
  geom_brain(atlas = dk, side = "lateral", hemi = "left")
p


set.seed(123)
dat <- ggplot_build(p)$data[[1]]
brain_region_list <- (unique(dat$group))
#brain_region_list <- brain_region_list[!is.na(brain_region_list)]
print(brain_region_list)

brain_region_point_list <- lapply(brain_region_list, function(x) {
  # print(x)
  pts <- dat[dat$group == x,]
  current_area <- st_area(pts$geometry)
  samples <- st_sample(pts$geometry, size = max(5, round(current_area/sampling_rate)), type= "regular", exact=F)
  # print(round(current_area/sampling_rate))
  tbl <- as_tibble(st_coordinates(samples)) %>% dplyr::mutate(label = x)
  tbl
})

brain_region_point <- bind_rows(brain_region_point_list) %>% dplyr::mutate(label = as.factor(label))
point_size_list <- sort(table(brain_region_point$label), decreasing = TRUE)
## Remove very small regions
# final_region_tbl <- brain_region_point %>% dplyr::filter(label %in% names(point_size_list)[1:15]) %>% dp

final_region_tbl <- brain_region_point %>%
  dplyr::filter(label %in% names(point_size_list)[1:15]) %>%
  dplyr::mutate(label = factor(label, levels = names(point_size_list)[1:15]))

final_region_tbl %>% ggplot(aes(x = X, y = Y, color = label)) + geom_point() + scale_color_discrete()
print(dim(brain_region_point))

# ---------------------------------------------------------------
library(scDesign3)
library(SingleCellExperiment)
library(tidyverse)
theme_set(theme_bw())
# example_sce <- readRDS((url("https://figshare.com/ndownloader/files/40581968")))
# saveRDS(~/R/graphBGM/example_sce, file = "example_sce.rds")

example_sce <- readRDS("~/R/graphBGM/example_sce.rds")
print(example_sce)

# Define a list of strings (a character vector)
my_list <- row.names(example_sce)
substring <- "ADT_"
matching_strings <- grep(substring, my_list, value = TRUE)
print(matching_strings)


cell_type_list <- table(colData(example_sce)$cell_type)
keep_gene <- c("CD4", "CD14", "CD19", "CD34", "CD3E", "CD8A")
keep_adt <- c("ADT_CD4", "ADT_CD14", "ADT_CD19", "ADT_CD34", "ADT_CD3", "ADT_CD8")
if (n_genes>800){
  keep_adt <- matching_strings
  print("Usa all ADT")
}
keep <- c(keep_gene, keep_adt)
idx <- which(rownames(example_sce) %in% keep)
idx <- c(1:(n_genes-length(keep_gene)),idx)
idx <- unique(c(1:(n_genes - length(keep_gene)), idx))
example_sce <- example_sce[idx,]
logcounts(example_sce) <- log1p(counts(example_sce))
rownames(example_sce) <- make.unique(rownames(example_sce))

example_data <- construct_data(
  sce = example_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "cell_type"
)

# Sys.setenv(BOOST_ROOT = "~/R/graphBGM/boost_1_82_0/")


# example_marginal <- fit_marginal(
#   data = example_data,
#   predictor = "gene",
#   mu_formula = "cell_type",
#   sigma_formula = "cell_type",
#   family_use = "nb",
#   n_cores = 1,
#   usebam = FALSE,
#   parallelization = "pbmcmapply"
# )

library(BiocParallel)
options(BiocParallel.packages = c("stats", "mgcv"))
library(stats)  # explicitly load 'stats'
library(mgcv)   # if using gam/bam models
param <- SnowParam(workers = 6, type = "SOCK")
register(param)
options(BiocParallel.loadPackages = TRUE)
example_marginal <- fit_marginal(
  data = example_data,
  predictor = "gene",
  mu_formula = "cell_type",
  sigma_formula = "cell_type",
  family_use = "nb",
  n_cores = 4,
  usebam = FALSE,
  parallelization = "bpmapply"
)
# example_marginal <- fit_marginal(
#   data = example_data,
#   predictor = "gene",
#   mu_formula = "cell_type",
#   sigma_formula = 0.001,
#   family_use = "nb",
#   n_cores = 4,
#   usebam = FALSE,
#   parallelization = "bpmapply"
# )

set.seed(123)
example_copula <- fit_copula(
  sce = example_sce,
  assay_use = "counts",
  marginal_list = example_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 1,
  input_data = example_data$dat,
  parallelization = "bpmapply"
)

### Set number of cells in cell type here
unique_combined <- example_data$dat %>% tidyr::expand(tidyr::nesting(cell_type, corr_group))
new_ct <- as.data.frame(unique_combined %>% uncount(as.numeric(point_size_list[1:15])))
example_para <- extract_para(
  sce = example_sce,
  marginal_list = example_marginal,
  n_cores = 1,
  family_use = "nb",
  new_covariate = new_ct,
  data = example_data$dat
)

set.seed(123)
print("Amplify the signal -- new here")
col_data <- new_ct %>% dplyr::bind_cols(final_region_tbl)
labels <- col_data$label
unique_labels <- unique(labels)
for (lbl in unique_labels) {
  idx <- which(labels == lbl)
  rand_val <- runif(1, min = 0, max = 1)
  example_para$mean_mat[idx, ] = 10* rand_val * example_para$mean_mat[idx, ] 
}


example_newcount <- simu_new(
  sce = example_sce,
  mean_mat = example_para$mean_mat,
  sigma_mat = example_para$sigma_mat,
  zero_mat = example_para$zero_mat,
  quantile_mat = NULL,
  copula_list = example_copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = example_data$dat,
  new_covariate = new_ct,
  important_feature = example_copula$important_feature,
  filtered_gene = example_data$filtered_gene
)


counts <- example_newcount
new_sce <- SingleCellExperiment(list(counts = counts), colData = col_data)

# saveRDS(new_sce, file = "simulated_mulimodalSTdata.rds")
g1 <- ggplot(col_data, aes(x=X, y=Y, color=label)) + geom_point()

g2 <- ggplot(col_data, aes(x=X, y=Y, color=cell_type)) + geom_point()

ggplot(col_data, aes(x=X, y=Y, color=corr_group)) + geom_point()

# Save files
df <- data.frame((as.matrix(counts)), row.names = counts@Dimnames[[1]])
features <- rownames(df)  
keep_rna <- setdiff(features, keep_adt)
exp_df <- t(df[keep_rna,])
adt_df <- t(df[keep_adt,])

n_cells <- dim(df)[2]
output_name <- paste0(data_name,'_nCells_', as.character(n_cells), '_nGenes_', as.character(n_genes))

out_dir <- ''
write.csv(exp_df, paste0(out_dir, 'Amplify_',output_name, '_RNA.csv'), row.names = FALSE)
write.csv(adt_df, paste0(out_dir, 'Amplify_',output_name, '_ADT.csv'), row.names = FALSE)
write.csv(col_data, paste0(out_dir, 'Amplify_',output_name, '_metadata.csv'), row.names = FALSE)
