args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
out_path <- paste(file, "_tmp/", sep = "")
out_path_reduction <- paste(out_path, "reduction/", sep = "")

data <- readRDS(file)

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}
if (!dir.exists(out_path_reduction)) {
  dir.create(out_path_reduction, recursive = TRUE)
}

dim_reds <- names(data@reductions)
for (dr in dim_reds) {
    rd <- data@reductions[[dr]]@cell.embeddings
    write.csv(rd, file = paste(out_path_reduction, dr, ".csv", sep = ""))
}

meta <- data@meta.data
write.csv(meta, file = paste(out_path, "meta.csv", sep = ""))

gene <- data@assays$RNA@meta.features
write.csv(gene, paste(out_path, "gene.csv", sep = ""))

counts <- data@assays$RNA@data
sparse_counts <- as(counts, "dgCMatrix")
Matrix::writeMM(sparse_counts, file = paste(out_path, "counts.mtx", sep = ""))