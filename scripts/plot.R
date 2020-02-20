suppressPackageStartupMessages(library(tidyverse))


filename_assess_assembly_all_meanQ_tsv <- "stats/assess_assembly_all_meanQ.tsv"
filename_assess_assembly_all_meanQ_pdf <- "assess_assembly_all_meanQ.pdf"

filename_assess_homopolymers_all_rel_len_tsv <- "stats/assess_homopolymers_all_rel_len.tsv"
filename_assess_homopolymers_all_rel_len_pdf <- "assess_homopolymers_all_rel_len.pdf"

filename_assess_homopolymers_all_correct_len_tsv <- "stats/assess_homopolymers_all_correct_len.tsv"
filename_assess_homopolymers_all_correct_len_pdf <- "assess_homopolymers_all_correct_len.pdf"


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
    theme_minimal() +
    theme(legend.position = "bottom") +
    ggtitle(paste("Mean Q-scores, reference", i_ref)) +
    ylab("") +
    xlab("") +
    theme(
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  filename_out <- paste0("plots/", i_ref, "_", filename_assess_assembly_all_meanQ_pdf)
  writeLines(paste("Saving plot to", filename_out))
  ggsave(p, filename = filename_out)
}


# plot assess_homopolymers correct vs lengths

df <- read_tsv(filename_assess_homopolymers_all_correct_len_tsv,
  col_names = c("assembly", "reference", "rlen", "A", "C", "G", "T", "A_n", "C_n", "G_n", "T_n", "AT_n", "AT", "GC_n", "GC"),
  col_types = "fffdddddddddddd"
) %>%
  select(-ends_with("_n")) %>%
  gather(key = "base", value = "frac_correct", -assembly, -reference, -rlen)

for (i_ref in unique(df$reference)) {
  p <- df %>%
    filter(reference == i_ref) %>%
    mutate(base = ordered(base, levels = c("A", "T", "C", "G", "AT", "GC"))) %>%
    group_by(assembly, reference, base, rlen) %>%
    ggplot(aes(x = rlen, y = frac_correct)) +
    geom_point(aes(color = assembly), na.rm=TRUE) +
    geom_line(aes(group = assembly, color = assembly), na.rm=TRUE) +
    facet_wrap(~base, ncol = 4, scales="free_x") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    ggtitle(paste("Fraction correct homopolymers, reference", i_ref)) +
    ylab("") + xlab("") +
    ylim(0, 1.0) +
    theme(
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90")
    ) +
    scale_color_discrete(guide = guide_legend(title="", ncol = 3))
  
  filename_out <- paste0("plots/", i_ref, "_", filename_assess_homopolymers_all_correct_len_pdf)
  writeLines(paste("Saving plot to", filename_out))
  ggsave(p, filename = filename_out)
}
