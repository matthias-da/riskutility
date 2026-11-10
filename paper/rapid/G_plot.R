## Figure for Experiment G: risk contribution of each quasi-identifier (drop-one, refitted),
## change in RAPID excess over the fixed-model null; points = mean, bars = min-max over releases.
suppressPackageStartupMessages(library(ggplot2))
here <- function(...) file.path("paper/rapid", ...)
res <- readRDS(here("results", "G_qi_contribution.rds"))
res <- res[res$dropped != "(none)", ]

agg <- do.call(rbind, lapply(split(res, list(res$target, res$dropped), drop = TRUE), function(d)
  data.frame(target = d$target[1], dropped = d$dropped[1], m = mean(d$d_excess),
             lo = min(d$d_excess), hi = max(d$d_excess))))
agg$target <- factor(agg$target, levels = c("marital.status", "income"),
                     labels = c("Sensitive: marital status", "Sensitive: income"))
lab <- function(v) {
  v <- sub("relationship+marital.status", "relationship + marital status", v, fixed = TRUE)
  v <- sub("education.num", "education", v, fixed = TRUE)
  gsub(".", " ", v, fixed = TRUE)
}
agg$group <- grepl("+", agg$dropped, fixed = TRUE)
agg <- agg[order(agg$target, agg$m), ]
agg$key <- factor(paste(agg$target, agg$dropped, sep = "___"),
                  levels = paste(agg$target, agg$dropped, sep = "___"))

p <- ggplot(agg, aes(x = m, y = key)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "grey55") +
  geom_point(aes(shape = group), size = 2, colour = "#1f4e79") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  scale_y_discrete(labels = function(k) lab(sub(".*___", "", k))) +
  facet_wrap(~ target, scales = "free", ncol = 2) +
  labs(x = "Drop in RAPID excess over the permutation null when the attacker does not know the QI",
       y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey95"))

out <- here("results", "plot_adult_qi_contribution.pdf")
ggsave(out, p, width = 7.5, height = 3.4)
cat("wrote", out, "\n")
print(agg[order(agg$target, -agg$m), c("target", "dropped", "m", "lo", "hi")], digits = 3,
      row.names = FALSE)
