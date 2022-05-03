suppressPackageStartupMessages(library(tidyverse))

filename_assess_assembly_all_scores_tsv <- "pomoxis/assess_assembly_all_scores.tsv"
filename_assess_assembly_all_meanQ_pdf <- "assess_assembly_all_meanQ.pdf"

# plot assess_assembly mean Q scores for each assembly and reference

df <- read_tsv(filename_assess_assembly_all_scores_tsv,
  col_names = c("assembly", "reference", "Qscore", "percErr"),
  col_types = "ffdc"
)

for (i_ref in unique(df$reference)) {
  df.plot <- df %>%
    filter(reference == i_ref) %>%
		#remove assemblies that were polished against other references
		filter(!(str_detect(assembly, "proovframe|homopolish") & !str_detect(assembly, i_ref)))

  p <- df.plot %>%
    ggplot(aes(x = reorder(assembly, Qscore, FUN = max), y = Qscore)) +
    geom_line(aes(group = reference), color = "grey70", size = 1) +
    geom_point(shape = 21, size = 3, fill = "deepskyblue3") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "Q-scores mean", subtitle = paste("Reference:", i_ref), caption = "from pomoxis/assess_assemblies output files *_summ.txt") +
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
