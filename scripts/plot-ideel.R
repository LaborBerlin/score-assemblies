suppressPackageStartupMessages(library(tidyverse))


# diamond for uniprot

dir_ideel_diamond <- "ideel/diamond/"
filename_ideel_pdf <- "ideel/ideel_stats.pdf"

filelist <- list.files(path = dir_ideel_diamond, pattern = ".*\\.tsv", recursive = FALSE, full.names = FALSE)
df_list <- vector("list", length(filelist))

for (i in seq_along(filelist)) {
  filename <- paste0(dir_ideel_diamond, filelist[[i]])
  df <- read_tsv(filename,
    col_names = c("qlen", "slen"),
    col_types = "ii"
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
  facet_wrap(~assembly, scales = "free_x") +
  theme_bw() +
  ggtitle(paste("ideel")) +
  theme(legend.position = "bottom") +
  ylab("") + xlab("qlen / slen") +
  xlim(0, 1.5) +
  scale_color_discrete(guide = FALSE) +
  scale_fill_discrete(guide = FALSE) +
  theme(strip.text.x = element_text(size = 6))

filename_out <- filename_ideel_pdf
writeLines(paste("Saving plot to", filename_out))
ggsave(p, filename = filename_out)


# diamond for references

dir_ideel_diamond_ref <- "ideel/diamond-ref/"

ref_list <- list.dirs(path = dir_ideel_diamond_ref, recursive = FALSE, full.names = FALSE)

for (r in seq_along(ref_list)) {


	filelist <- list.files(path = paste0(dir_ideel_diamond_ref,"/", ref_list[[r]]), pattern = ".*\\.tsv", recursive = FALSE, full.names = FALSE)
	df_list <- vector("list", length(filelist))

	for (i in seq_along(filelist)) {
		filename <- paste0(dir_ideel_diamond_ref, "/", ref_list[[r]], "/", filelist[[i]])
		df <- read_tsv(filename,
			col_names = c("qlen", "slen"),
			col_types = "ii"
		) %>%
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
		ggtitle(paste("ideel")) +
		theme(legend.position = "bottom") +
		ylab("") + xlab("qlen / slen") +
		xlim(0, 1.5) +
		scale_color_discrete(guide = FALSE) +
		scale_fill_discrete(guide = FALSE) +
		theme(strip.text.x = element_text(size = 4))

	filename_pdf <- paste0("ideel/", ref_list[[r]], "_ideel_stats.pdf")
	writeLines(paste("Saving plot to", filename_pdf))
	ggsave(p, filename = filename_pdf)

	p <- df.all %>%
		filter(slen == qlen) %>%
		group_by(assembly) %>%
			summarize(n = n()) %>%
		ungroup() %>%
		ggplot(aes(x=reorder(assembly,n), y=n)) +
		geom_point() +
		geom_text(aes(label=n), vjust = 1.5, size = 2.5) +
		coord_flip() +
		xlab("") +
		labs(title=ref_list[[r]], subtitle = "number of proteins with qlen = slen") +
		theme_bw()

	filename_pdf <- paste0("ideel/", ref_list[[r]], "_samelen.pdf")
	writeLines(paste("Saving plot to", filename_pdf))
	ggsave(p, filename = filename_pdf)


}
