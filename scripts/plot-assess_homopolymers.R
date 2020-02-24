suppressPackageStartupMessages(library(tidyverse))

filename_assess_homopolymers_all_correct_len_tsv <- "pomoxis/assess_homopolymers_all_correct_len.tsv"
filename_assess_homopolymers_all_correct_len_pdf <- "assess_homopolymers_all_correct_len.pdf"

# plot assess_homopolymers correct vs lengths

df <- read_tsv(filename_assess_homopolymers_all_correct_len_tsv,
  col_names = c("assembly", "reference", "rlen", "A", "C", "G", "T", "A_n", "C_n", "G_n", "T_n", "AT_n", "AT", "GC_n", "GC"),
  col_types = "fffdddddddddddd"
) %>%
  select(-ends_with("_n")) %>%
  gather(key = "base", value = "frac_correct", -assembly, -reference, -rlen)
df
for (i_ref in unique(df$reference)) {
  p <- df %>%
    filter(reference == i_ref) %>%
    mutate(base = ordered(base, levels = c("A", "T", "C", "G", "AT", "GC"))) %>%
    ggplot(aes(x = rlen, y = frac_correct)) +
    geom_point(aes(color = assembly), na.rm = TRUE) +
    geom_line(aes(group = assembly, color = assembly), na.rm = TRUE) +
    facet_wrap(~base, ncol = 4, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "Fraction correct homopolymers", subtitle = paste("Reference:", i_ref), caption = "from pomoxis/assess_homopolymers output files *_correct_len.tsv") +
    ylab("") +
    xlab("") +
    ylim(0, 1.0) +
    theme(
      strip.background = element_rect(fill = "grey90")
    ) +
    scale_color_discrete(guide = guide_legend(title = "", ncol = 3))

  filename_out <- paste0("pomoxis/", i_ref, "_", filename_assess_homopolymers_all_correct_len_pdf)
  writeLines(paste("Saving plot to", filename_out))
  ggsave(p, filename = filename_out)
}

