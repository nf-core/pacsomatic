#!/usr/bin/env Rscript

library(DSS)
library(bsseq)
library(data.table)

# Get input files and parameters from Nextflow
tumor_file <- "${tumor_bed}"
normal_file <- "${normal_bed}"
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
output_file <- paste0(prefix, ".dmr.tsv")
threads <- as.integer("${task.cpus}")

# Parse additional arguments if provided
args <- "${args}"
p_threshold <- 0.001
delta <- 0.25
minCG <- 5
minlen <- 100
smoothing_span <- 500

# Parse args string if provided (format: "p_threshold=0.001 delta=0.25 minCG=5 minlen=100 smoothing_span=500")
if (args != "") {
    arg_list <- strsplit(args, " ")[[1]]
    for (arg in arg_list) {
        if (grepl("=", arg)) {
            parts <- strsplit(arg, "=")[[1]]
            if (parts[1] == "p_threshold") p_threshold <- as.numeric(parts[2])
            if (parts[1] == "delta") delta <- as.numeric(parts[2])
            if (parts[1] == "minCG") minCG <- as.integer(parts[2])
            if (parts[1] == "minlen") minlen <- as.integer(parts[2])
            if (parts[1] == "smoothing_span") smoothing_span <- as.integer(parts[2])
        }
    }
}

# Log parameters
cat("Running DSS DMR analysis with parameters:\\n")
cat("  Tumor file:", tumor_file, "\\n")
cat("  Normal file:", normal_file, "\\n")
cat("  Output file:", output_file, "\\n")
cat("  Threads:", threads, "\\n")
cat("  p_threshold:", p_threshold, "\\n")
cat("  delta:", delta, "\\n")
cat("  minCG:", minCG, "\\n")
cat("  minlen:", minlen, "\\n")
cat("  smoothing_span:", smoothing_span, "\\n\\n")

# Prepare input files
# Handle gzipped files
if (grepl("\\\\.gz\$", tumor_file)) {
    tumor_tmp <- tempfile(pattern = "tumor_", fileext = ".bed")
    system2("gunzip", args = c("-c", tumor_file), stdout = tumor_tmp)
    tumor_tmp1 <- tempfile(pattern = "tumor1_", fileext = ".bed")
    system2("grep", args = c("-v '^#'", tumor_tmp), stdout = tumor_tmp1)
    tumor_file <- tumor_tmp1
}

if (grepl("\\\\.gz\$", normal_file)) {
    normal_tmp <- tempfile(pattern = "normal_", fileext = ".bed")
    system2("gunzip", args = c("-c", normal_file), stdout = normal_tmp)
    normal_tmp1 <- tempfile(pattern = "normal1_", fileext = ".bed")
    system2("grep", args = c("-v '^#'", normal_tmp), stdout = normal_tmp1)
    normal_file <- normal_tmp1
}

# Extract required columns (chr, pos, N, X) if file has more columns
# Assuming columns 1,2,6,7 based on the original script
cat("Reading tumor file...\\n")
tumor_raw <- fread(tumor_file, header = FALSE, sep = "\\t")
if (ncol(tumor_raw) >= 7) {
    tumor <- tumor_raw[, .(V1, V2, V6, V7)]
} else if (ncol(tumor_raw) == 4) {
    tumor <- tumor_raw
} else {
    stop("Unexpected number of columns in tumor file. Expected 4 or more columns.")
}
colnames(tumor) <- c("chr", "pos", "N", "X")

cat("Reading normal file...\\n")
normal_raw <- fread(normal_file, header = FALSE, sep = "\\t")
if (ncol(normal_raw) >= 7) {
    normal <- normal_raw[, .(V1, V2, V6, V7)]
} else if (ncol(normal_raw) == 4) {
    normal <- normal_raw
} else {
    stop("Unexpected number of columns in normal file. Expected 4 or more columns.")
}
colnames(normal) <- c("chr", "pos", "N", "X")

rm(tumor_raw, normal_raw)

cat("Creating BSseq object...\\n")
bs_obj <- makeBSseqData(
    list(
        tumor,
        normal
    ),
    sampleNames = c("tumor", "normal")
)

rm(list = c("tumor", "normal"))
gc()

cat("Performing DML test with smoothing...\\n")
dmlTest.sm <- DMLtest(
    bs_obj,
    group1 = c("tumor"),
    group2 = c("normal"),
    smoothing = TRUE,
    smoothing.span = smoothing_span,
    ncores = threads
)

rm(list = c("bs_obj"))
gc()

cat("Calling DMRs...\\n")
dmrs <- callDMR(
    dmlTest.sm,
    p.threshold = p_threshold,
    delta = delta,
    minCG = minCG,
    minlen = minlen
)

# Check if DMRs were found
if (is.null(dmrs) || nrow(dmrs) == 0) {
    cat("Warning: No DMRs found with current thresholds.\\n")
    # Create empty output file with header
    empty_df <- data.frame(
        chr = character(),
        start = integer(),
        end = integer(),
        length = integer(),
        nCG = integer(),
        meanMethy1 = numeric(),
        meanMethy2 = numeric(),
        diff.Methy = numeric(),
        areaStat = numeric()
    )
    fwrite(
        empty_df,
        output_file,
        sep = "\\t",
        quote = FALSE,
        na = "NA",
        col.names = TRUE
    )
} else {
    cat("Found", nrow(dmrs), "DMRs\\n")
    fwrite(
        dmrs,
        output_file,
        sep = "\\t",
        quote = FALSE,
        na = "NA",
        col.names = TRUE,
        nThread = threads
    )
}

# Clean up temporary files if created
if (exists("tumor_tmp") && file.exists(tumor_tmp)) {
    unlink(tumor_tmp)
}
if (exists("normal_tmp") && file.exists(normal_tmp)) {
    unlink(normal_tmp)
}

# Print session info
cat("\\nSession Info:\\n")
sessionInfo()

# Write versions file
writeLines(
    c(
        '"${task.process}":',
        paste('    dss:', as.character(packageVersion('DSS'))),
        paste('    bsseq:', as.character(packageVersion('bsseq'))),
        paste('    data.table:', as.character(packageVersion('data.table')))
    ),
'versions.yml')

cat("\\nDSS DMR analysis completed successfully!\\n")
