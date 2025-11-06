#!/usr/bin/env Rscript

library(BiocManager)
ref_genome <- "${ref_genome}"
if (!require(ref_genome, quietly = TRUE)) {
  local_lib <- file.path(getwd(), ".library")
  if (!dir.exists(local_lib)) {
    dir.create(local_lib, recursive = TRUE)
  }
  .libPaths(c(local_lib, .libPaths()))

  BiocManager::install(ref_genome, update = FALSE, ask = FALSE)
}

library(MutationalPatterns)
library(ref_genome, character.only = TRUE)
library(cowplot)
library(ggplot2)

set.seed(42)

vc_file <- "${vcf}"
sample_name <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
max_delta <- as.numeric("${max_delta}")

vcf_grl <- read_vcfs_as_granges(vc_file, sample_name, ref_genome)

muts <- mutations_from_vcf(vcf_grl[[1]])

# types <- mut_type(vcf_grl[[1]])
# context <- mut_context(vcf_grl[[1]], ref_genome)
# type_context <- type_context(vcf_grl[[1]], ref_genome)

type_occurences <- mut_type_occurrences(vcf_grl[[1]], ref_genome)
# Write to file
write.table(type_occurences, file = paste0(sample_name, ".type_occurences.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

mut_mat <- mut_matrix(vcf_list = vcf_grl, ref_genome = ref_genome)

signatures = get_known_signatures()

fit_res_bs <- fit_to_signatures_bootstrapped(mut_mat, signatures, max_delta = max_delta, n_boots = 50, method = "strict")
fit_res <- fit_to_signatures_strict(mut_mat, signatures, max_delta = max_delta)
rownames(fit_res\$fit_res\$reconstructed) <- rownames(mut_mat)

# Write bootstrapped mutsig
write.table(
    fit_res_bs,
    file = paste0(sample_name, ".mut_sigs_bootstrapped.tsv"),
    sep = "\\t",
    quote = FALSE,
    row.names = TRUE, col.names = NA
    )

# Write mut sig
write.table(
    fit_res\$fit_res\$contribution[order(fit_res\$fit_res\$contribution, decreasing = TRUE),,drop=FALSE],
    file = paste0(sample_name, ".mut_sigs.tsv"),
    sep = "\\t",
    quote = FALSE,
    row.names = TRUE, col.names = NA
    )

# Write reconstructed sig
write.table(
    fit_res\$fit_res\$reconstructed,
    file = paste0(sample_name, ".reconstructed_sigs.tsv"),
    sep = "\\t",
    quote = FALSE,
    row.names = TRUE, col.names = NA
    )

# Plot all
cont_plot <- plot_contribution(contribution = fit_res\$fit_res\$contribution, signatures)
profile_96_plot <- plot_96_profile(mut_mat)
cosine_sim_plot <- plot_original_vs_reconstructed(mut_mat, fit_res\$fit_res\$reconstructed)

combined_plot <- plot_grid(
    cosine_sim_plot +
        geom_text(aes(y = cos_sim, label = paste0(round(cos_sim, 3) * 100, "%")), vjust = 0) +
        scale_y_continuous(breaks = c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0)),
    cont_plot +
        geom_text(aes(label = ifelse(Contribution > 0, round(Contribution/sum(Contribution), 2), "")), position = position_fill(vjust = 0.5)),
    profile_96_plot,
    ncol = 3,
    rel_widths = c(0.1, 0.2, 0.7),
    align = "free"
    )

ggsave(paste0(sample_name, ".mutation_profile.pdf"), combined_plot, width=16, height=8, useDingbats=FALSE)

# versions file
writeLines(
    c(
        '"${task.process}":',
        paste('    mutationalpattern:', as.character(packageVersion('MutationalPatterns'))),
        paste('    cowplot:', as.character(packageVersion('cowplot'))),
        paste('    ggplot2:', as.character(packageVersion('ggplot2')))
    ),
'versions.yml')
