suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidytext))

dir <- ifelse(!is.na(snakemake@params[['out_dir']]), snakemake@params[['out_dir']], '.')

filename_nucdiff_stats_tsv <- paste0(dir, "/nucdiff/all_stats.tsv")
filename_suffix_nucdiff_stats_tsv <- "_nucdiff_stats.pdf"

# plot nucdiff stats for each reference

df <- read_tsv(filename_nucdiff_stats_tsv,
  col_names = c("assembly", "reference", "measure", "value"),
  col_types = "fffd"
) %>%
  mutate(measure = fct_relevel(measure, "Substitutions"))

for (i_ref in unique(df$reference)) {
  df.plot <- df %>%
    filter(reference == i_ref)
  
  p <- df.plot %>%
    filter(reference == i_ref) %>%
    group_by(measure) %>% 
      mutate(assembly = fct_reorder(assembly, -value)) %>%
    ungroup() %>%
    ggplot(aes(y = assembly, x = value)) +
    #geom_line(aes(group = reference), color = "grey70", size = 1) +
    geom_point(shape = 21, size = 2, fill = "deepskyblue3") +
    facet_wrap(~measure, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "nucdiff", subtitle = paste("Reference:", i_ref)) +
    ylab("") +
    xlab("") +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text = element_text(size = rel(0.75))
    )
  filename_out <- paste0(dir, "/nucdiff/", i_ref, filename_suffix_nucdiff_stats_tsv)
  writeLines(paste("Saving plot to", filename_out))
  ggsave(p, filename = filename_out)
}

