suppressPackageStartupMessages(library(tidyverse))

dir <- ifelse(!is.na(snakemake@params[['out_dir']]), snakemake@params[['out_dir']], '.')

# diamond for uniprot -------------------------------------------------------------------------------------
dir_ideel_diamond <- paste0(dir, "/ideel/diamond/")
filename_ideel_hist_pdf <- paste0(dir, "/ideel/ideel_uniprot_histograms.pdf")
filename_ideel_box_pdf <- paste0(dir, "/ideel/ideel_uniprot_boxplots.pdf")

filelist <- list.files(path = dir_ideel_diamond, pattern = ".*\\.tsv", recursive = FALSE, full.names = FALSE)
df_list <- vector("list", length(filelist))

for (i in seq_along(filelist)) {
  filename <- paste0(dir_ideel_diamond, filelist[[i]])
	df <- read_tsv(filename,
									 col_names = c("qseqid", "sseqid", "qlen", "slen"),
									 col_types = "ccii"
    ) %>%
    mutate(assembly = filelist[[i]]) %>%
    mutate(assembly = str_remove(assembly, "\\.tsv")) %>%
    mutate(SCov = qlen / slen)

  df_list[[i]] <- df
}
df.all <- dplyr::bind_rows(df_list) %>% mutate(assembly = factor(assembly))

p <- ggplot(df.all, aes(x = SCov)) +
  geom_histogram(aes(color=assembly, fill=assembly), alpha = 0.8, position = "identity", binwidth = 0.01, size=0.2) +
  #geom_density(aes(fill = assembly), alpha = 0.4, color = "grey40") +
  facet_wrap(~assembly) +
  theme_bw() +
  ggtitle(paste("ideel with Uniprot Sprot")) +
  ylab("") + xlab("qlen / slen") +
  xlim(0.5, 1.5) +
  scale_color_discrete(guide = FALSE) +
  scale_fill_discrete(guide = FALSE) +
  theme(strip.text.x = element_text(size = 6)) +
  theme(axis.text = element_text(size = rel(0.75)))

filename_out <- filename_ideel_hist_pdf
writeLines(paste("Saving plot to", filename_out))
ggsave(p, filename = filename_out)


p.box <- df.all %>%
  mutate(assembly = fct_reorder(assembly, SCov, mean)) %>%
  ggplot(aes(y = assembly, x = SCov)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_histogram(aes(color=assembly, fill=assembly), alpha = 0.8, position = "identity", binwidth = 0.01, size=0.2) +
  #geom_density(aes(fill = assembly), alpha = 0.4, color = "grey40") +
  #facet_wrap(~assembly) +
  theme_bw() +
  ggtitle(paste("ideel with Uniprot Sprot")) +
  ylab("") + xlab("qlen / slen") +
  xlim(0.5, 1.5) +
  theme(axis.text = element_text(size = rel(0.75)))

filename_out <- filename_ideel_box_pdf
writeLines(paste("Saving plot to", filename_out))
ggsave(p.box, filename = filename_out)


# diamond for references --------------------------------------------------------------------------------
dir_ideel_diamond_ref <- paste0(dir, "/ideel/diamond-ref/")

ref_list <- list.dirs(path = dir_ideel_diamond_ref, recursive = FALSE, full.names = FALSE)

for (r in seq_along(ref_list)) {
  filelist <- list.files(path = paste0(dir_ideel_diamond_ref,"/", ref_list[[r]]), pattern = ".*\\.tsv", recursive = FALSE, full.names = FALSE)
  df_list <- vector("list", length(filelist))

  for (i in seq_along(filelist)) {
    filename <- paste0(dir_ideel_diamond_ref, "/", ref_list[[r]], "/", filelist[[i]])
    df <- read_tsv(filename, col_names = c("qseqid", "sseqid", "qlen", "slen"), col_types = "ccii") %>%
      mutate(assembly = filelist[[i]]) %>%
      mutate(assembly = str_remove(assembly, "\\.tsv")) %>%
      mutate(assembly = str_remove(assembly, paste0("_",ref_list[[r]]))) %>%
      mutate(SCov = qlen / slen)

    df_list[[i]] <- df
  }
  df.all <- dplyr::bind_rows(df_list) %>% mutate(assembly = factor(assembly))

	p <- ggplot(df.all, aes(x = SCov)) +
		geom_histogram(aes(color=assembly, fill=assembly), alpha = 0.8, position = "identity", binwidth = 0.01, size=0.2) +
		#geom_density(aes(fill = assembly), alpha = 0.4, color = "grey40") +
		facet_wrap(~assembly, scales = "free_x") +
		theme_bw() +
		labs(title = "ideel", subtitle = paste("Reference:", ref_list[[r]])) +
		ylab("") + xlab("qlen / slen") +
		xlim(0.5, 1.5) +
		scale_color_discrete(guide = FALSE) +
		scale_fill_discrete(guide = FALSE) +
		theme(strip.text.x = element_text(size = 6)) +
		theme(axis.text = element_text(size = rel(0.75)))

	filename_pdf <- paste0(dir, "/ideel/", ref_list[[r]], "_ideel_histograms.pdf")
	writeLines(paste("Saving plot to", filename_pdf))
	ggsave(p, filename = filename_pdf)

	p.box <- df.all %>%
	  mutate(assembly = fct_reorder(assembly, SCov, mean)) %>%
	  ggplot(aes(y = assembly, x = SCov)) +
	  geom_boxplot(outlier.shape = NA) +
	  theme_bw() +
	  labs(title = "ideel", subtitle = paste("Reference:", ref_list[[r]])) +
	  ylab("") + xlab("qlen / slen") +
	  xlim(0.5, 1.5) +
	  theme(axis.text.y = element_text(size = rel(0.75)))
	
	filename_pdf <- paste0(dir, "/ideel/", ref_list[[r]], "_ideel_boxplots.pdf")
	writeLines(paste("Saving plot to", filename_pdf))
	ggsave(p.box, filename = filename_pdf)
	
}
