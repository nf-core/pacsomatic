#!/usr/bin/env Rscript

library(annotatr)
library(GenomicRanges)
library(data.table)

# Get input files and parameters from Nextflow
dmr_file <- "${dmr}"
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
threads <- as.integer("${task.cpus}")

# Parse additional arguments if provided
args <- "${args}"
genome <- "hg38"
ignore_strand <- TRUE

# Parse args string if provided (format: "genome=hg38 ignore_strand=TRUE")
if (args != "") {
    arg_list <- strsplit(args, " ")[[1]]
    for (arg in arg_list) {
        if (grepl("=", arg)) {
            parts <- strsplit(arg, "=")[[1]]
            if (parts[1] == "genome") genome <- parts[2]
            if (parts[1] == "ignore_strand") ignore_strand <- as.logical(parts[2])
        }
    }
}

# Log parameters
cat("Running annotatr DMR annotation with parameters:\\n")
cat("  DMR file:", dmr_file, "\\n")
cat("  Prefix:", prefix, "\\n")
cat("  Genome:", genome, "\\n")
cat("  Ignore strand:", ignore_strand, "\\n")
cat("  Threads:", threads, "\\n\\n")

# Read DMR file
cat("Reading DMR file...\\n")
dmrs <- fread(dmr_file, nThread = threads)

# Check if DMRs file is empty or has no data
if (nrow(dmrs) == 0) {
    cat("Warning: No DMRs found in input file.\\n")

    # Create empty summary file
    dm_annsum <- data.table(
        Annotation_Type = character(),
        Count = integer()
    )
    fwrite(
        dm_annsum,
        paste0(prefix, "_dmr_annotation_summary.tsv.gz"),
        sep = "\\t",
        quote = FALSE,
        na = "NA",
        compress = "gzip",
        nThread = threads
    )

    cat("Created empty annotation summary file.\\n")

} else {
    cat("Found", nrow(dmrs), "DMRs\\n")

    # Convert to GRanges
    cat("Converting to GRanges...\\n")
    dmrs_gr <- GRanges(dmrs)

    # Define annotations to build
    ann_to_build <- c(
        paste0(genome, "_genes_promoters"),
        paste0(genome, "_genes_1to5kb"),
        paste0(genome, "_genes_5UTRs"),
        paste0(genome, "_genes_exons"),
        paste0(genome, "_genes_introns"),
        paste0(genome, "_genes_3UTRs")
    )

    cat("Building annotations for genome:", genome, "\\n")
    cat("Annotations:", paste(ann_to_build, collapse = ", "), "\\n")

    # Build annotations
    built_annotations <- build_annotations(
        genome = genome,
        annotations = ann_to_build
    )

    # Annotate DMRs
    cat("Annotating DMRs...\\n")
    dmrs_annotated <- annotate_regions(
        dmrs_gr,
        built_annotations,
        ignore.strand = ignore_strand
    )

    # Summarize annotation type
    cat("Summarizing annotations...\\n")
    dm_annsum <- summarize_annotations(
        annotated_regions = dmrs_annotated,
        quiet = TRUE
    )

    colnames(dm_annsum) <- c("Annotation_Type", "Count")

    # Write annotation summary
    cat("Writing annotation summary...\\n")
    fwrite(
        dm_annsum,
        paste0(prefix, "_dmr_annotation_summary.tsv.gz"),
        sep = "\\t",
        quote = FALSE,
        na = "NA",
        append = FALSE,
        compress = "gzip",
        nThread = threads
    )

    # Split by annotation type and write each type to a file
    cat("Writing annotation-specific files...\\n")
    for (ann_type in ann_to_build) {
        # Extract DMRs for this annotation type
        # Check if annotation exists in the data
        if (ann_type %in% dmrs_annotated\$annot\$type) {
            ann_type_dmrs <- dmrs_annotated[dmrs_annotated\$annot\$type == ann_type]

            # Convert to data.table
            ann_type_dt <- as.data.table(ann_type_dmrs)

            # Only write if there are DMRs for this annotation
            if (nrow(ann_type_dt) > 0) {
                fwrite(
                    ann_type_dt,
                    paste0(prefix, "_", ann_type, "_dmrs.tsv.gz"),
                    sep = "\\t",
                    quote = FALSE,
                    na = "NA",
                    append = FALSE,
                    compress = "gzip",
                    nThread = threads
                )
                cat("  Wrote", nrow(ann_type_dt), "DMRs for", ann_type, "\\n")
            } else {
                cat("  No DMRs found for", ann_type, "\\n")
            }
        } else {
            cat("  Annotation type", ann_type, "not found in results\\n")
        }
    }

    cat("\\nAnnotation complete!\\n")
}

# Print session info
cat("\\nSession Info:\\n")
sessionInfo()

# Write versions file
writeLines(
    c(
        '"${task.process}":',
        paste('    annotatr:', as.character(packageVersion('annotatr'))),
        paste('    genomicranges:', as.character(packageVersion('GenomicRanges'))),
        paste('    data.table:', as.character(packageVersion('data.table')))
    ),
'versions.yml')

cat("\\nDMR annotation completed successfully!\\n")
