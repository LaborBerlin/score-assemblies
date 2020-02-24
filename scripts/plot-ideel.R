suppressPackageStartupMessages(library(tidyverse))

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
  scale_fill_discrete(guide = FALSE)

filename_out <- filename_ideel_pdf
writeLines(paste("Saving plot to", filename_out))
ggsave(p, filename = filename_out)

