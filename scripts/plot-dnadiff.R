suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidytext))

filename_dnadiff_stats_tsv <- "dnadiff/all_stats.tsv"
filename_suffix_dnadiff_stats_tsv <- "_dnadiff_stats.pdf"

# plot dnadiff stats for each reference

df <- read_tsv(filename_dnadiff_stats_tsv,
  col_names = c("assembly", "reference", "measure", "value", "value2"),
  col_types = "fffdd"
)

for (i_ref in unique(df$reference)) {
  df.plot <- df %>%
    filter(reference == i_ref)
  
  p <- df.plot %>%
    group_by(measure) %>% 
      mutate(assembly = fct_reorder(assembly, value)) %>%
    ungroup() %>%
    ggplot(aes(y = assembly, x = value)) +
    geom_point(shape = 21, size = 2, fill = "deepskyblue3") +
    facet_wrap(~ measure, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "dnadiff", subtitle = paste("Reference:", i_ref)) +
    ylab("") +
    xlab("") +
    theme(axis.text.y = element_text(size = rel(0.75)))

  filename_out <- paste0("dnadiff/", i_ref, filename_suffix_dnadiff_stats_tsv)
  writeLines(paste("Saving plot to", filename_out))
  ggsave(p, filename = filename_out)
}

