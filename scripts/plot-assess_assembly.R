suppressPackageStartupMessages(library(tidyverse))

filename_assess_assembly_all_meanQ_tsv <- "pomoxis/assess_assembly_all_meanQ.tsv"
filename_assess_assembly_all_meanQ_pdf <- "assess_assembly_all_meanQ.pdf"

# plot assess_assembly mean Q scores for each assembly and reference

df <- read_tsv(filename_assess_assembly_all_meanQ_tsv,
  col_names = c("assembly", "reference", "Qscore"),
  col_types = "ffd"
)

for (i_ref in unique(df$reference)) {
  p <- df %>%
    filter(reference == i_ref) %>%
    ggplot(aes(x = reorder(assembly, Qscore, FUN = max), y = Qscore)) +
    geom_line(aes(group = reference), color = "grey70", size = 1) +
    geom_point(shape = 21, size = 3, fill = "deepskyblue3") +
    theme_bw() +
    theme(legend.position = "bottom") +
    ggtitle(paste("Mean Q-scores, reference", i_ref)) +
    ylab("") +
    xlab("") +
    theme(
      strip.background = element_rect(fill = "grey90"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  filename_out <- paste0("pomoxis/", i_ref, "_", filename_assess_assembly_all_meanQ_pdf)
  writeLines(paste("Saving plot to", filename_out))
  ggsave(p, filename = filename_out)
}

