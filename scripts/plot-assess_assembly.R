suppressPackageStartupMessages(library(tidyverse))

dir <- ifelse(!is.na(snakemake@params[['out_dir']]), snakemake@params[['out_dir']], '.')

filename_assess_assembly_all_scores_tsv <- paste0(dir, "/pomoxis/assess_assembly_all_scores.tsv")
filename_assess_assembly_all_meanQ_pdf <- "assess_assembly_all_meanQ.pdf"

# plot assess_assembly mean Q scores for each assembly and reference

df <- read_tsv(filename_assess_assembly_all_scores_tsv,
  col_names = c("assembly", "reference", "Qscore", "percErr"),
  col_types = "ffdc"
)

for (i_ref in unique(df$reference)) {
  df.plot <- df %>%
    filter(reference == i_ref)

  p <- df.plot %>%
    ggplot(aes(y = reorder(assembly, Qscore, FUN = max), x = Qscore)) +
    geom_point(shape = 21, size = 2, fill = "deepskyblue3") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "pomoxis / assess-assembly", subtitle = paste("Reference:", i_ref)) +
    ylab("") +
    xlab("Q-score ") +
    theme(axis.text = element_text(size = rel(0.75)))
  
  filename_out <- paste0(dir, "/pomoxis/", i_ref, "_", filename_assess_assembly_all_meanQ_pdf)
  writeLines(paste("Saving plot to", filename_out))
  ggsave(p, filename = filename_out)
}
