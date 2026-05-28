exc <- read.delim("SEAAD_DLPCF/merged.tpm.Exc.txt", check.names = FALSE)
inh <- read.delim("SEAAD_DLPCF/merged.tpm.Inh.txt", check.names = FALSE)

gene_exc <- exc[[1]]
gene_inh <- inh[[1]]

exc_values <- exc[, -1]
inh_values <- inh[, -1]

exc_mean <- rowMeans(exc_values, na.rm = TRUE)
inh_mean <- rowMeans(inh_values, na.rm = TRUE)

exc_df <- data.frame(gene = gene_exc, excitatory_mean = exc_mean)
inh_df <- data.frame(gene = gene_inh, inhibitory_mean = inh_mean)

plot_df <- merge(exc_df, inh_df, by = "gene")

write.csv(plot_df, "exc_vs_inh_gene_mean_TPM.csv", row.names = FALSE)

png("exc_vs_inh_scatter_rawTPM.png", width = 1800, height = 1500, res = 300)
plot(plot_df$excitatory_mean,
     plot_df$inhibitory_mean,
     pch = 16,
     cex = 0.4,
     xlab = "Mean TPM in excitatory neurons",
     ylab = "Mean TPM in inhibitory neurons",
     main = "Gene Expression Comparison: Excitatory vs Inhibitory Neurons")
abline(0, 1, lty = 2)
dev.off()

png("exc_vs_inh_scatter_logTPM.png", width = 1800, height = 1500, res = 300)
plot(log10(plot_df$excitatory_mean + 1),
     log10(plot_df$inhibitory_mean + 1),
     pch = 16,
     cex = 0.4,
     xlab = "log10(mean TPM in excitatory neurons + 1)",
     ylab = "log10(mean TPM in inhibitory neurons + 1)",
     main = "Gene Expression Comparison: Excitatory vs Inhibitory Neurons")
abline(0, 1, lty = 2)
dev.off()
