suppressPackageStartupMessages(library(tidyverse))

filename_busco_stats_tsv <- "busco/all_stats.tsv"
filename_busco_pdf <- "busco/busco_stats.pdf"

# plot busco

df <- read_tsv(filename_busco_stats_tsv,
  col_names = c("assembly", "Complete", "Fragmented", "Missing", "n"),
  col_types = "fdddd"
) %>%
  gather(key = "measure", value = "percent", -assembly, -n)

p <- df %>%
  ggplot(aes(y = reorder(assembly, percent, max), x = percent)) +
  geom_text(data = filter(df, measure == "Complete"), aes(color = measure, label = percent), size = 2.5, hjust=1.4) +
  geom_point(size = 2, aes(color = measure)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "BUSCO") +
  xlab("Percent") +
  ylab("") +
  xlim(0, 100) +
  theme(axis.text.y = element_text(size = rel(0.75))) +
  scale_color_discrete(guide = guide_legend(title = ""))

filename_out <- filename_busco_pdf
writeLines(paste("Saving plot to", filename_out))
ggsave(p, filename = filename_out)

