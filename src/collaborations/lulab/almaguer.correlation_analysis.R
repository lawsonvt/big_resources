library(ggplot2)
library(cowplot)
library(cocor)

# read in data
set1 <- read.csv("~/Documents/projects/lulab/berenice_almaguer/distance_km_v_ctrl_r.txt")
set2 <- read.csv("~/Documents/projects/lulab/berenice_almaguer/distance_km_v_cko_r.txt")



# correlation analyses

cor1 <- cor.test(set1$ctrl_r, set1$distance_km, method = "pearson")
cor2 <- cor.test(set2$cko_r, set2$distance_km, method = "pearson")

# Display the correlations
cat("Cell Type 1 correlation: r =", round(cor1$estimate, 3), 
    ", p =", round(cor1$p.value, 4), "\n")
cat("Cell Type 2 correlation: r =", round(cor2$estimate, 3), 
    ", p =", round(cor2$p.value, 4), "\n")

# plot em
p1 <- ggplot(set1, aes(x = distance_km, y = ctrl_r)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", color = "darkblue", se = TRUE) +
  theme_bw() +
  annotate("text", x = Inf, y = Inf, 
           label = paste("r =", round(cor1$estimate, 3),
                         "\np =", round(cor1$p.value, 3)),
           hjust = 1.1, vjust = 1.5, size = 5)


p2 <- ggplot(set2, aes(x = distance_km, y = cko_r)) +
  geom_point(color = "coral", size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  theme_bw() +
  annotate("text", x = Inf, y = Inf, 
           label = paste("r =", round(cor2$estimate, 3),
                         "\np =", round(cor2$p.value, 3)),
           hjust = 1.1, vjust = 1.5, size = 5)

plot_grid(p1, p2)
ggsave("~/Documents/projects/lulab/berenice_almaguer/correlation_plots.png", width=9, height=4)

cor_comparison <- data.frame(
  cell_type = c("Type 1", "Type 2"),
  correlation = c(cor1$estimate, cor2$estimate),
  ci_lower = c(cor1$conf.int[1], cor2$conf.int[1]),
  ci_upper = c(cor1$conf.int[2], cor2$conf.int[2])
)

ggplot(cor_comparison, aes(x = cell_type, y = correlation, fill = cell_type)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  scale_fill_manual(values = c("Type 1" = "steelblue", "Type 2" = "coral")) +
  labs(title = "Correlation Comparison with 95% Confidence Intervals",
       x = "Cell Type",
       y = "Pearson Correlation (r)") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")


ggplot(cor_comparison, aes(x = correlation, y = cell_type, color = cell_type)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, size = 1) +
  scale_color_manual(values = c("Type 1" = "steelblue", "Type 2" = "coral")) +
  labs(title = "Correlation Coefficients with 95% CIs",
       x = "Pearson Correlation (r)",
       y = "Cell Type") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")

# For independent groups (different mice for each cell type)
result <- cocor.indep.groups(
  r1.jk = cor(set1$ctrl_r, set1$distance_km),  # correlation for set1
  r2.hm = cor(set2$cko_r, set2$distance_km),  # correlation for set2
  n1 = nrow(set1),                               # sample size set1
  n2 = nrow(set2),                               # sample size set2
  alternative = "two.sided",
  test = "all"  # runs multiple comparison tests
)

# Extract values
r1 <- cor(set1$ctrl_r, set1$distance_km)
r2 <- cor(set2$cko_r, set2$distance_km)
n1 <- nrow(set1)
n2 <- nrow(set2)

# Fisher z-transformation
z1 <- 0.5 * log((1 + r1) / (1 - r1))
z2 <- 0.5 * log((1 + r2) / (1 - r2))

# 1. Visualize the transformation from r to z space
r_values <- seq(-0.99, 0.99, by = 0.01)
z_values <- 0.5 * log((1 + r_values) / (1 - r_values))

transformation_df <- data.frame(r = r_values, z = z_values)

ggplot(transformation_df, aes(x = r, y = z)) +
  geom_line(color = "gray30", size = 1) +
  geom_point(aes(x = r1, y = z1), color = "steelblue", size = 5) +
  geom_point(aes(x = r2, y = z2), color = "coral", size = 5) +
  geom_segment(aes(x = r1, y = z1, xend = r1, yend = -4), 
               linetype = "dashed", color = "steelblue", alpha = 0.5) +
  geom_segment(aes(x = r2, y = z2, xend = r2, yend = -4), 
               linetype = "dashed", color = "coral", alpha = 0.5) +
  annotate("text", x = r1, y = z1 + 0.3, label = "Type 1", color = "steelblue", size = 4) +
  annotate("text", x = r2, y = z2 + 0.3, label = "Type 2", color = "coral", size = 4) +
  labs(title = "Fisher z-Transformation",
       subtitle = "Converting correlations to z-space for comparison",
       x = "Pearson r (correlation)",
       y = "Fisher z") +
  theme_minimal()

# 2. Compare correlations in both r-space and z-space
comparison_data <- data.frame(
  cell_type = rep(c("Type 1", "Type 2"), 2),
  value = c(r1, r2, z1, z2),
  space = rep(c("r-space", "z-space"), each = 2)
)

ggplot(comparison_data, aes(x = cell_type, y = value, fill = cell_type)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  facet_wrap(~space, scales = "free_y") +
  scale_fill_manual(values = c("Type 1" = "steelblue", "Type 2" = "coral")) +
  labs(title = "Correlations in r-space vs z-space",
       subtitle = "Fisher transformation normalizes the distribution",
       x = "Cell Type",
       y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")


p_value <- result@fisher1925$p.value
z_stat <- (z1 - z2) / sqrt(1/(n1-3) + 1/(n2-3))

# Create a standard normal distribution
x_norm <- seq(-4, 4, by = 0.01)
y_norm <- dnorm(x_norm)

norm_df <- data.frame(x = x_norm, y = y_norm)

# Critical values for alpha = 0.05
critical_value <- qnorm(0.975)

ggplot(norm_df, aes(x = x, y = y)) +
  geom_line(size = 1) +
  geom_area(data = subset(norm_df, abs(x) > critical_value), 
            aes(x = x, y = y), fill = "red", alpha = 0.3) +
  geom_vline(xintercept = z_stat, color = "darkblue", size = 1.5, linetype = "dashed") +
  geom_vline(xintercept = c(-critical_value, critical_value), 
             color = "red", linetype = "dotted") +
  annotate("text", x = z_stat, y = max(y_norm) * 0.8, 
           label = paste("Z =", round(z_stat, 3)), 
           color = "darkblue", size = 5, hjust = -0.1) +
  annotate("text", x = 0, y = max(y_norm) * 0.5, 
           label = paste("p =", round(p_value, 4)), 
           color = "darkblue", size = 5) +
  labs(title = "Fisher z-Test: Test Statistic Distribution",
       subtitle = "Standard normal distribution with critical regions (α = 0.05)",
       x = "z-statistic",
       y = "Density") +
  theme_minimal()

# 4. Effect size visualization - difference between correlations
diff_r <- abs(r1 - r2)
diff_z <- abs(z1 - z2)

effect_df <- data.frame(
  metric = c("Difference in r", "Difference in z", "Z-statistic"),
  value = c(diff_r, diff_z, abs(z_stat))
)

ggplot(effect_df, aes(x = metric, y = value, fill = metric)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Effect Size Metrics",
       subtitle = paste("p-value =", round(p_value, 4)),
       x = "",
       y = "Value") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))

# 5. Create a comprehensive summary figure
library(gridExtra)
library(grid)

# Text summary
summary_text <- paste0(
  "Fisher z-Test Results\n",
  "─────────────────────\n",
  "Type 1: r = ", round(r1, 3), " (n = ", n1, ")\n",
  "Type 2: r = ", round(r2, 3), " (n = ", n2, ")\n",
  "\n",
  "Fisher z-values:\n",
  "Type 1: z = ", round(z1, 3), "\n",
  "Type 2: z = ", round(z2, 3), "\n",
  "\n",
  "Test statistic: Z = ", round(z_stat, 3), "\n",
  "P-value: ", round(p_value, 4), "\n",
  "\n",
  ifelse(p_value < 0.05, 
         "Conclusion: Correlations differ significantly", 
         "Conclusion: No significant difference")
)

text_grob <- textGrob(summary_text, 
                      gp = gpar(fontfamily = "mono", fontsize = 10),
                      just = "left")

grid.arrange(text_grob)

# new appproach, fit a linear model

# Combine datasets with group indicator
combined_data <- rbind(
  data.frame(cell_count = set1$ctrl_r, 
             distance = set1$distance_km, 
             cell_type = "Type1"),
  data.frame(cell_count = set2$cko_r, 
             distance = set2$distance_km, 
             cell_type = "Type2")
)

# Fit model with interaction
model <- lm(distance ~ cell_count * cell_type, data = combined_data)
summary(model)

# The interaction term tests if slopes differ
# Look at: cell_count:cell_typeType2

# Visualize
ggplot(combined_data, aes(x = cell_count, y = distance, color = cell_type)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("Type1" = "steelblue", "Type2" = "coral")) +
  labs(title = "Testing for Different Slopes via Interaction",
       subtitle = paste("Interaction p-value:", 
                        round(summary(model)$coefficients["cell_count:cell_typeType2", "Pr(>|t|)"], 4))) +
  theme_minimal()
ggsave("~/Documents/projects/lulab/berenice_almaguer/linear_model_fit.png", width=7, height=6)

# permutation test

# Function to calculate difference in correlations
calc_cor_diff <- function(data1, data2) {
  r1 <- cor(data1$cell_count, data1$distance)
  r2 <- cor(data2$cell_count, data2$distance)
  return(abs(r1 - r2))
}

# Observed difference
observed_diff <- calc_cor_diff(combined_data[combined_data$cell_type == "Type1",], 
                               combined_data[combined_data$cell_type == "Type2",])

# Permutation test
set.seed(123)
n_permutations <- 10000


n1 <- nrow(set1)
n_total <- nrow(combined_data)

# Run permutations
perm_diffs <- replicate(n_permutations, {
  # Randomly shuffle and split
  shuffled_idx <- sample(1:n_total)
  perm_set1 <- combined_data[shuffled_idx[1:n1], ]
  perm_set2 <- combined_data[shuffled_idx[(n1+1):n_total], ]
  
  calc_cor_diff(perm_set1, perm_set2)
})

# Calculate p-value
p_value_perm <- mean(perm_diffs >= observed_diff)

# Visualize permutation distribution
hist(perm_diffs, breaks = 50, col = "lightblue", 
     main = "Permutation Test Distribution",
     xlab = "Difference in Correlations",
     xlim = c(0, max(c(perm_diffs, observed_diff))))
abline(v = observed_diff, col = "red", lwd = 2)
text(observed_diff, max(table(cut(perm_diffs, 50))), 
     paste("Observed\np =", round(p_value_perm, 4)), 
     pos = 4, col = "red")

