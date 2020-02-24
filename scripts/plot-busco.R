suppressPackageStartupMessages(library(tidyverse))

filename_busco_stats_tsv <- "busco/all_stats.tsv"
filename_busco_pdf <- "busco/busco_stats.pdf"

# plot busco

df <- read_tsv(filename_busco_stats_tsv,
               col_names = c("assembly", "Complete", "Fragmented", "Missing", "n"),
               col_types = "fdddd") %>% 
  gather(key = "measure", value = "percent", -assembly, -n)

p <- df %>%
    ggplot(aes(x = reorder(assembly, percent, max), y = percent)) +
    geom_line(aes(group = measure, color = measure), size = 1) +
    geom_point(size = 2, aes(color = measure)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    ggtitle(paste("BUSCO")) +
    ylab("") +
    xlab("") +
    theme(
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_color_discrete(guide = guide_legend(title=""))

filename_out <- filename_busco_pdf
writeLines(paste("Saving plot to", filename_out))
ggsave(p, filename = filename_out)
