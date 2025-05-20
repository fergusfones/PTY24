seeds <- 1:10

results_list <- list()
cor_list <- list()
pval_list <- list()

mouseCTD <- ewceData::ctd()

#Loop over each seed

for (i in seeds) {
  set.seed(i)  # Set the seed
  module_df <- data.frame(
    gene_id = names(netwk$colors),
    colors = labels2colors(netwk$colors)
  )

  # Get Module Eigengenes per cluster
  MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

  # Reorder modules so similar modules are next to each other
  MEs0 <- orderMEs(MEs0)
  module_order = names(MEs0) %>% gsub("ME","", .)

  # Add treatment names
  MEs0$treatment = row.names(MEs0)


  # Adding other variables
  MEs0 <- cbind(MEs0, meta_df_CA1_neun)

  ### Remove grey module (junk)
  MEs0 <- MEs0[,-which(colnames(MEs0)=="MEgrey")]

  traits <- meta_df_CA1_neun %>%
    select(- c("SegmentDisplayName", "Histology.no.", "population", "grossRegion"))

  #encoding traits
  # do as factor
  traits$Group <- as.factor(traits$Group)

  traits_numeric <- ifelse(traits == "WW", 0 ,1)

  # module trait correlation
  # adjust rows in MEs0 to just show modules
  mod_cols <- grep("^ME", colnames(MEs0), value = TRUE)
  moduleTraitCor = stats::cor(MEs0[, mod_cols],traits_numeric , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs0))

  # Store the correlation and p-values
  cor_list[[paste0("seed_", i)]] <- moduleTraitCor
  pval_list[[paste0("seed_", i)]] <- moduleTraitPvalue


  # Split modules
  gene_lists <- split(module_df$gene_id, module_df$colors)

    # Loop over all modules
    for (mod_name in names(gene_lists)) {
      genes <- gene_lists[[mod_name]]

      # Run EWCE
      results <- EWCE::bootstrap_enrichment_test(
        sct_data = mouseCTD,
        sctSpecies = "mouse",
        genelistSpecies = "mouse",
        hits = genes,
        reps = 1000,
        annotLevel = 1
      )

      # Save results with metadata
      res_df <- results$results
      res_df$seed <- i
      res_df$module <- mod_name
      res_df$n_genes <- length(genes)
      results_list[[paste0("seed", i, "_", mod_name)]] <- res_df
  }

  # Save intermediate results after each seed
  saveRDS(results_list, file = paste0("results_list_seed_", i, ".rds"))
  saveRDS(cor_list, file = paste0("cor_list_seed_", i, ".rds"))
  saveRDS(pval_list, file = paste0("pval_list_seed_", i, ".rds"))


}



# Compare the results (example: correlation matrix, overlap, etc.)
# This part depends on the format of your 'results' – e.g., enrichment p-values or scores.



# multi - plot of all seeds

ewce_combined <- bind_rows(results_list)

ggplot(ewce_combined, aes(x = CellType, y = sd_from_mean, fill = as.factor(seed))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ seed) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot the heatmap


# Convert correlation list to long data frame
cor_long <- map_dfr(names(cor_list), function(seed) {
  cor_df <- as.data.frame(cor_list[[seed]])
  cor_df$Module <- rownames(cor_df)
  cor_df %>%
    pivot_longer(-Module, names_to = "Trait", values_to = "Correlation") %>%
    mutate(Seed = seed)
})

# Do the same for p-values
pval_long <- map_dfr(names(pval_list), function(seed) {
  pval_df <- as.data.frame(pval_list[[seed]])
  pval_df$Module <- rownames(pval_df)
  pval_df %>%
    pivot_longer(-Module, names_to = "Trait", values_to = "Pvalue") %>%
    mutate(Seed = seed)
})

# Merge them for plotting
cor_pval_long <- left_join(cor_long, pval_long, by = c("Seed", "Module", "Trait"))


ggplot(cor_pval_long, aes(x = Module, y = Correlation, fill = Trait)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Seed) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module-Trait Correlations Across Seeds")



module_sizes <- map_dfr(names(results_list), function(name) {
  df <- results_list[[name]]
  data.frame(
    seed = df$seed[1],
    module = df$module[1],
    n_genes = df$n_genes[1]
  )
})

ggplot(module_sizes, aes(x = module, y = n_genes, fill = as.factor(seed))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + facet_wrap(~seed)


