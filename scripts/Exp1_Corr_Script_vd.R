# ============================================================================
# EXPERIMENT I (correlational) 
# ============================================================================
library(tidyverse)
library(readr)
library(psych)
library(rsq)
library(ggtext)
library(BayesFactor)
library(mediation)
library(broom)
library(gt)
library(grid)
library(conflicted)

conflict_prefer("select",  "dplyr")
conflict_prefer("filter",  "dplyr")
conflict_prefer("recode",  "dplyr")
conflict_prefer("mediate", "mediation")

# ============================================================================
# 1. Load data
# ============================================================================
df_main <- read_csv("data/data_experiment_1_corr_Main_vf  - dat.csv")
df_ctrl <- read_csv("data/data_experiment_1_corr_wControl_vf - dat.csv")

# ============================================================================
# 2. Rename variables
# ============================================================================
df_main <- df_main %>%
  rename(
    happy_post          = happy_post_1,
    happy_pre           = happy_pre_1,
    perceived_comp_post = understanding_post_1,
    perceived_comp_pre  = understanding_pre_1,
    desire_learn        = avg_learn
  ) %>%
  mutate(domain = factor(domain, levels = c("art", "phil")))   # Art History = reference

df_ctrl <- df_ctrl %>%
  mutate(
    across(-c(RecordedDate, ResponseId, domain, Condition), as.numeric),
    domain    = factor(domain, levels = c("art", "phil")),
    Condition = factor(Condition, levels = c("AI", "No"))       # AI feedback = reference
  )

# ============================================================================
# 3. PCA on the five AI-tone ratings -> Perceived sycophantic tone
# ============================================================================
items5 <- c("empathetic_1", "polite_1", "helpful_1", "blunt_1", "confrontational_1")
tone5  <- df_main %>%
  select(all_of(items5)) %>%
  na.omit() %>%
  mutate(across(everything(), as.numeric))

# Sampling adequacy
KMO(tone5)
cortest.bartlett(cor(tone5), n = nrow(tone5))

# Parallel analysis 
n <- nrow(tone5); p <- ncol(tone5)
obs_eigen <- eigen(cor(tone5))$values

set.seed(2026)
n_iter  <- 1000
eig_sim <- matrix(NA_real_, n_iter, p)   # simulated random-normal null
eig_res <- matrix(NA_real_, n_iter, p)   # permutation null (real marginals)
for (i in seq_len(n_iter)) {
  eig_sim[i, ] <- eigen(cor(matrix(rnorm(n * p), n, p)))$values
  eig_res[i, ] <- eigen(cor(apply(tone5, 2, sample)))$values
}
pa_table <- data.frame(
  component = seq_len(p),
  observed  = obs_eigen,
  sim_p95   = apply(eig_sim, 2, quantile, probs = .95),
  res_p95   = apply(eig_res, 2, quantile, probs = .95)
)
print(pa_table, digits = 3)

# Velicer's MAP
map_vals <- VSS(tone5, n = 4, rotate = "none", fm = "pc", plot = FALSE)$map
cat("MAP minimum at", which.min(map_vals), "component(s)\n")
print(round(map_vals, 3))

# Loadings and variance 
pca_result <- principal(tone5, nfactors = 1, rotate = "none")
print(pca_result$loadings, cutoff = 0, digits = 3)
pca_result$Vaccounted

# Component scores
pc_fit     <- prcomp(tone5, center = TRUE, scale. = TRUE)
pc1_scores <- pc_fit$x[, 1]
if (cor(pc1_scores, tone5$empathetic_1) < 0) pc1_scores <- -pc1_scores

df_main$sycophancy_score <- NA_real_
cc <- complete.cases(df_main[, items5])
df_main$sycophancy_score[cc] <- pc1_scores

# ============================================================================
# 4. Z-score all continuous variables used in modelling
# ============================================================================
z_vars <- c("sycophancy_score",
            "happy_pre", "happy_post",
            "perceived_comp_pre", "perceived_comp_post",
            "confi_pre", "confi_post",
            "accuracy_pre", "accuracy_post",
            "desire_learn")

df_main <- df_main %>%
  mutate(across(all_of(z_vars), ~ as.numeric(scale(.))))

# Z-score control dataset
df_ctrl <- df_ctrl %>%
  mutate(across(any_of(c("accuracy_pre", "accuracy_post")), ~ as.numeric(scale(.))))

# ============================================================================
# 5. PRIMARY MODELS 
# ============================================================================

## S1.1 Mood
model_S1_1a <- lm(happy_post ~ sycophancy_score + happy_pre, data = df_main)
summary(model_S1_1a)
model_S1_1b <- lm(happy_post ~ sycophancy_score + happy_pre +
                    accuracy_pre + accuracy_post + domain, data = df_main)
summary(model_S1_1b)

## S1.2 Perceived competence
model_S1_2a <- lm(perceived_comp_post ~ sycophancy_score + perceived_comp_pre, data = df_main)
summary(model_S1_2a)
model_S1_2b <- lm(perceived_comp_post ~ sycophancy_score + perceived_comp_pre +
                    accuracy_pre + accuracy_post + domain, data = df_main)
summary(model_S1_2b)

## S1.3 Confidence 
model_S1_3a <- lm(confi_post ~ sycophancy_score + confi_pre, data = df_main)
summary(model_S1_3a)
model_S1_3b <- lm(confi_post ~ sycophancy_score + confi_pre +
                    accuracy_pre + accuracy_post + domain, data = df_main)
summary(model_S1_3b)

## S1.4 Desire to learn
model_S1_4a <- lm(desire_learn ~ sycophancy_score, data = df_main)
summary(model_S1_4a)
model_S1_4b <- lm(desire_learn ~ sycophancy_score +
                    accuracy_pre + accuracy_post + domain, data = df_main)
summary(model_S1_4b)

## S1.5 Actual competence 
model_S1_5a <- lm(accuracy_post ~ sycophancy_score, data = df_main)
summary(model_S1_5a)
model_S1_5b <- lm(accuracy_post ~ sycophancy_score + accuracy_pre + domain, data = df_main)
summary(model_S1_5b)

## S1.6 Final exam performance: AI feedback vs. no-feedback control (df_ctrl)
model_S1_6a <- lm(accuracy_post ~ Condition, data = df_ctrl)
summary(model_S1_6a)
model_S1_6b <- lm(accuracy_post ~ Condition + accuracy_pre + domain, data = df_ctrl)
summary(model_S1_6b)

# ============================================================================
# 6. BAYES FACTORS  (BF01)
# ============================================================================
df_bf <- as.data.frame(df_main)
df_bf$domain <- factor(df_bf$domain)
df_bf <- df_bf[complete.cases(df_bf[, c("accuracy_post", "sycophancy_score",
                                        "accuracy_pre", "domain")]), ]

## Actual competence
set.seed(2026)
bf_syco     <- lmBF(accuracy_post ~ sycophancy_score, data = df_bf)
bf01_nocont <- 1 / bf_syco
bf01_nocont

## Actual competence — with controls 
set.seed(2026)
bf_full    <- lmBF(accuracy_post ~ sycophancy_score + accuracy_pre + domain, data = df_bf)
bf_reduced <- lmBF(accuracy_post ~ accuracy_pre + domain,                    data = df_bf)
bf01_cont  <- bf_reduced / bf_full
bf01_cont

# ============================================================================
# 7. MEDIATION MODELS  
# ============================================================================

## S2.1 Perceived competence
# a) 
med_m_S2_1a <- lm(happy_post ~ sycophancy_score + happy_pre, data = df_main)
out_m_S2_1a <- lm(perceived_comp_post ~ sycophancy_score + happy_post + happy_pre +
                    perceived_comp_pre, data = df_main)
set.seed(2026)
med_S2_1a <- mediate(med_m_S2_1a, out_m_S2_1a,
                     treat = "sycophancy_score", mediator = "happy_post",
                     boot = TRUE, sims = 1000)
summary(med_S2_1a)

# b) 
med_m_S2_1b <- lm(happy_post ~ sycophancy_score + happy_pre +
                    accuracy_pre + accuracy_post + domain, data = df_main)
out_m_S2_1b <- lm(perceived_comp_post ~ sycophancy_score + happy_post + happy_pre +
                    perceived_comp_pre + accuracy_pre + accuracy_post + domain, data = df_main)
set.seed(2026)
med_S2_1b <- mediate(med_m_S2_1b, out_m_S2_1b,
                     treat = "sycophancy_score", mediator = "happy_post",
                     boot = TRUE, sims = 1000)
summary(med_S2_1b)

## S2.2 Confidence
# a) 
med_m_S2_2a <- lm(happy_post ~ sycophancy_score + happy_pre, data = df_main)
out_m_S2_2a <- lm(confi_post ~ sycophancy_score + happy_post + happy_pre +
                    confi_pre, data = df_main)
set.seed(2026)
med_S2_2a <- mediate(med_m_S2_2a, out_m_S2_2a,
                     treat = "sycophancy_score", mediator = "happy_post",
                     boot = TRUE, sims = 1000)
summary(med_S2_2a)

# b) 
med_m_S2_2b <- lm(happy_post ~ sycophancy_score + happy_pre +
                    accuracy_pre + accuracy_post + domain, data = df_main)
out_m_S2_2b <- lm(confi_post ~ sycophancy_score + happy_post + happy_pre +
                    confi_pre + accuracy_pre + accuracy_post + domain, data = df_main)
set.seed(2026)
med_S2_2b <- mediate(med_m_S2_2b, out_m_S2_2b,
                     treat = "sycophancy_score", mediator = "happy_post",
                     boot = TRUE, sims = 1000)
summary(med_S2_2b)

## S2.3 Desire to learn
# a) 
med_m_S2_3a <- lm(happy_post ~ sycophancy_score + happy_pre, data = df_main)
out_m_S2_3a <- lm(desire_learn ~ sycophancy_score + happy_post + happy_pre, data = df_main)
set.seed(2026)
med_S2_3a <- mediate(med_m_S2_3a, out_m_S2_3a,
                     treat = "sycophancy_score", mediator = "happy_post",
                     boot = TRUE, sims = 1000)
summary(med_S2_3a)

# b) 
med_m_S2_3b <- lm(happy_post ~ sycophancy_score + happy_pre +
                    accuracy_pre + accuracy_post + domain, data = df_main)
out_m_S2_3b <- lm(desire_learn ~ sycophancy_score + happy_post + happy_pre +
                    accuracy_pre + accuracy_post + domain, data = df_main)
set.seed(2026)
med_S2_3b <- mediate(med_m_S2_3b, out_m_S2_3b,
                     treat = "sycophancy_score", mediator = "happy_post",
                     boot = TRUE, sims = 1000)
summary(med_S2_3b)



# ============================================================================
# FIGURES 
# ============================================================================

tone_labels_map <- c(empathetic_1 = "Empathy", polite_1 = "Politeness",
                     helpful_1 = "Helpfulness", blunt_1 = "Bluntness",
                     confrontational_1 = "Confrontation")
tone_levels <- c("Empathy", "Politeness", "Helpfulness", "Bluntness", "Confrontation")
tone_colors <- c(Empathy = "#FF0000", Politeness = "#FF8C00", Helpfulness = "#228B22",
                 Bluntness = "#87CEEB", Confrontation = "#00008B")

## ---- Fig 2a: tone-ratings violin plot (raw 0-100 scale) ----
# NOTE: recompute the raw ratings here because df_main tone items are not z-scored,
# but for safety pull the raw values before any scaling if this block is re-run.
df_tone_long <- df_main %>%
  pivot_longer(cols = names(tone_labels_map), names_to = "Tone", values_to = "Rating") %>%
  mutate(Tone = factor(tone_labels_map[Tone], levels = tone_levels))

ggplot(df_tone_long, aes(x = Tone, y = Rating, fill = Tone)) +
  geom_violin(trim = FALSE, alpha = 0.4, scale = "width") +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.6, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5, fill = "white") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8) +
  scale_fill_manual(values = tone_colors) +
  labs(x = "", y = "Ratings") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 30),
        legend.position = "none")

## ---- Fig 2b: PCA loadings (PC1) ----
loadings_df <- data.frame(
  Variable = rownames(unclass(pca_result$loadings)),
  Loading  = as.numeric(unclass(pca_result$loadings)[, 1])
) %>%
  mutate(Variable = recode(Variable, !!!tone_labels_map),
         Variable = factor(Variable, levels = tone_levels))

ggplot(loadings_df, aes(x = Loading, y = Variable)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(Loading, 2)), hjust = 0.5, size = 4) +
  coord_flip(xlim = c(-0.4, 0.9)) +
  labs(x = "Factor Loadings", y = NULL) +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(size = 30, angle = 35, hjust = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

## ---- Fig 3a-e: scatter plots (z-scored axes) ----
scatter_theme <- function() {
  theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title.x = element_markdown(size = 30, lineheight = 0.1),
          axis.title.y = element_markdown(size = 30, lineheight = 0.1))
}
x_lab <- paste0("<b style='font-weight:900;'>Perceived sycophantic tone</b>",
                "<br><span style='font-size:20pt;'>(PCA score)</span>")
z_lim <- c(-3, 3)

# Fig 3a - Mood
ggplot(df_main, aes(x = sycophancy_score, y = happy_post)) +
  geom_jitter(width = 0.05, height = 0.05, alpha = 0.35, size = 0.6, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "orange", fill = "orange", alpha = 0.2) +
  coord_cartesian(ylim = z_lim) +
  labs(x = x_lab, y = "<b style='font-weight:900;'>Mood</b><br><span style='font-size:20pt;'>(post AI-interaction, z-scored)</span>") +
  scatter_theme()

# Fig 3b - Confidence
ggplot(df_main, aes(x = sycophancy_score, y = confi_post)) +
  geom_jitter(width = 0.05, height = 0.05, alpha = 0.35, size = 0.6, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "lightblue", fill = "lightblue", alpha = 0.2) +
  coord_cartesian(ylim = z_lim) +
  labs(x = x_lab, y = "<b style='font-weight:900;'>Confidence</b><br><span style='font-size:20pt;'>(post AI-interaction, z-scored)</span>") +
  scatter_theme()

# Fig 3c - Desire to learn
ggplot(df_main, aes(x = sycophancy_score, y = desire_learn)) +
  geom_jitter(width = 0.05, height = 0.05, alpha = 0.35, size = 0.6, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "forestgreen", fill = "forestgreen", alpha = 0.2) +
  coord_cartesian(ylim = z_lim) +
  labs(x = x_lab, y = "<b style='font-weight:900;'>Desire to learn</b><br><span style='font-size:20pt;'>(post AI-interaction, z-scored)</span>") +
  scatter_theme()

# Fig 3d - Perceived competence
ggplot(df_main, aes(x = sycophancy_score, y = perceived_comp_post)) +
  geom_jitter(width = 0.05, height = 0.05, alpha = 0.35, size = 0.6, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "blue", alpha = 0.2) +
  coord_cartesian(ylim = z_lim) +
  labs(x = x_lab, y = "<b style='font-weight:900;'>Perceived competence</b><br><span style='font-size:20pt;'>(post AI-interaction, z-scored)</span>") +
  scatter_theme()

# Fig 3e - Actual competence
ggplot(df_main, aes(x = sycophancy_score, y = accuracy_post)) +
  geom_jitter(width = 0.05, height = 0.05, alpha = 0.35, size = 0.6, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.2) +
  coord_cartesian(ylim = z_lim) +
  labs(x = x_lab,
       y = paste0("<b style='font-weight:900;'>Competence</b><br>",
                  "<span style='font-size:20pt;'>(final exam score, z-scored)</span>")) +
  scatter_theme()



