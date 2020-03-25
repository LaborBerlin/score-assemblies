suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidytext))

filename_nucdiff_stats_tsv <- "nucdiff/all_stats.tsv"
filename_suffix_nucdiff_stats_tsv <- "_nucdiff_stats.pdf"

# plot nucdiff stats for each reference

df <- read_tsv(filename_nucdiff_stats_tsv,
  col_names = c("assembly", "reference", "measure", "value"),
  col_types = "fffd"
)

for (i_ref in unique(df$reference)) {
  p <- df %>%
    filter(reference == i_ref) %>%
    mutate(assembly = reorder_within(assembly, value, measure, FUN = min)) %>%
    ggplot(aes(x = assembly, y = value)) +
    geom_line(aes(group = reference), color = "grey70", size = 1) +
    geom_point(shape = 21, size = 3, fill = "deepskyblue3") +
    scale_x_reordered() +
    facet_wrap(~measure, scales = "free") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "nucdiff", subtitle = paste("Reference:", i_ref), caption = "from nucdiff output files *.report") +
    ylab("") +
    xlab("") +
    theme(
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  filename_out <- paste0("nucdiff/", i_ref, filename_suffix_nucdiff_stats_tsv)
  writeLines(paste("Saving plot to", filename_out))
  ggsave(p, filename = filename_out)
}

