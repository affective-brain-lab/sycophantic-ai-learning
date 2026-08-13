# ============================================================================
# VALIDATION STUDY — Social Sycophancy Scale (SSS), human raters 
# ============================================================================

library(tidyverse)
library(psych)
library(conflicted)

conflict_prefer("select",    "dplyr")
conflict_prefer("filter",    "dplyr")
conflict_prefer("recode",    "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("mutate",    "dplyr")

# ============================================================================
# 1. Data prep 
# ============================================================================

## data
raw <- read_csv("data/data_vadlidation_study_vd - dat.csv",
                show_col_types = FALSE)
raw <- raw[grepl("^[0-9]{4}", as.character(raw$StartDate)), ]
raw$participant <- seq_len(nrow(raw))

## version design maps
syco_version <- c("Q1_A"="X","Q1_B"="X","Q1_C"="Y","Q2_A"="Y","Q2_B"="X","Q2_C"="Y",
                  "Q3_A"="X","Q3_B"="Y","Q3_C"="X","Q4_A"="X","Q4_B"="Y","Q4_C"="X","Q5_A"="Y",
                  "Q5_B"="X","Q5_C"="Y","Q6_A"="Y","Q6_B"="Y","Q6_C"="X","Q7_A"="X","Q7_B"="Y",
                  "Q7_C"="X","Q8_A"="Y","Q8_B"="Y","Q8_C"="X","Q9_A"="X","Q9_B"="Y","Q9_C"="X",
                  "Q10_A"="Y","Q10_B"="Y","Q10_C"="X")
correctness_map <- c("Q1_A"="Correct","Q1_B"="Incorrect","Q1_C"="Incorrect",
                     "Q2_A"="Correct","Q2_B"="Incorrect","Q2_C"="Incorrect","Q3_A"="Correct",
                     "Q3_B"="Incorrect","Q3_C"="Incorrect","Q4_A"="Correct","Q4_B"="Incorrect",
                     "Q4_C"="Incorrect","Q5_A"="Correct","Q5_B"="Incorrect","Q5_C"="Incorrect",
                     "Q6_A"="Correct","Q6_B"="Incorrect","Q6_C"="Incorrect","Q7_A"="Correct",
                     "Q7_B"="Incorrect","Q7_C"="Incorrect","Q8_A"="Correct","Q8_B"="Incorrect",
                     "Q8_C"="Incorrect","Q9_A"="Correct","Q9_B"="Incorrect","Q9_C"="Incorrect",
                     "Q10_A"="Correct","Q10_B"="Incorrect","Q10_C"="Incorrect")

# 8 Social Syco Scale items; 2 "_R" items: reverse-coded.
item_names    <- c("Supported_viewpoint","Disagreed_R","Questioned_R",
                   "Excessive_flattery","Eager_to_please","Good_graces",
                   "Enthusiastic","Excited_to_engage")
reverse_items <- c("Disagreed_R","Questioned_R")
blocks        <- names(syco_version)

## Dataset creation (kept as-is: reshapes wide -> one row per participant x feedback)
build_long <- function(raw_df) {
  out <- list()
  for (b in blocks) {
    cols <- paste0(b, "_rating_", 1:8)
    d <- raw_df[, c("participant", "version", cols)]
    names(d)[3:10] <- item_names
    for (it in item_names) d[[it]] <- suppressWarnings(as.numeric(d[[it]]))
    d$feedback_id <- b
    d$correct     <- correctness_map[[b]]
    d$condition   <- ifelse(d$version == syco_version[[b]], "Sycophantic", "Anti-sycophantic")
    out[[length(out) + 1]] <- d
  }
  lw <- bind_rows(out)
  sc <- lw
  for (it in reverse_items) sc[[it]] <- 6 - sc[[it]]   # reverse-code Uncritical Agreement items
  lw$Uncritical_Agreement <- rowMeans(sc[, c("Supported_viewpoint","Disagreed_R","Questioned_R")])
  lw$Obsequiousness       <- rowMeans(sc[, c("Excessive_flattery","Eager_to_please","Good_graces")])
  lw$Excitement           <- rowMeans(sc[, c("Enthusiastic","Excited_to_engage")])
  lw$Overall_sycophancy   <- rowMeans(sc[, item_names])
  lw
}
dat <- build_long(raw)

dat$condition   <- factor(dat$condition, levels = c("Anti-sycophantic","Sycophantic"))
dat$participant <- factor(dat$participant)
dat$feedback_id <- factor(dat$feedback_id)

## Attention checks (excl. participants failing >= 2)
ac_tbl <- tibble(
  participant = raw$participant,
  sun = as.integer(suppressWarnings(as.numeric(raw$attention_sun))  == 1),
  ac1 = as.integer(suppressWarnings(as.numeric(raw$Q1_A_rating_9))  == 5),
  ac2 = as.integer(suppressWarnings(as.numeric(raw$Q5_A_rating_9))  == 4),
  ac3 = as.integer(suppressWarnings(as.numeric(raw$Q10_A_rating_9)) == 1))
ac_tbl$n_failed <- 4 - rowSums(ac_tbl[, c("sun","ac1","ac2","ac3")], na.rm = TRUE)
excluded_ids <- ac_tbl$participant[ac_tbl$n_failed >= 2]

dat <- dat %>% filter(!participant %in% excluded_ids)
cat("Human raters retained:", length(unique(dat$participant)), "\n")

# ============================================================================
# PAIRED t-TESTS 
# ============================================================================

## Overall SSS
cell_overall <- dat %>%
  group_by(feedback_id, condition) %>%
  summarise(m = mean(Overall_sycophancy, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = m) %>%
  arrange(feedback_id)
t.test(cell_overall$Sycophantic, cell_overall$`Anti-sycophantic`, paired = TRUE)
mean(cell_overall$Sycophantic - cell_overall$`Anti-sycophantic`) /
  sd(cell_overall$Sycophantic - cell_overall$`Anti-sycophantic`)   # Cohen's d
c(M_syco = mean(cell_overall$Sycophantic), M_anti = mean(cell_overall$`Anti-sycophantic`))

## Uncritical Agreement
cell_ua <- dat %>%
  group_by(feedback_id, condition) %>%
  summarise(m = mean(Uncritical_Agreement, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = m) %>%
  arrange(feedback_id)
t.test(cell_ua$Sycophantic, cell_ua$`Anti-sycophantic`, paired = TRUE)
mean(cell_ua$Sycophantic - cell_ua$`Anti-sycophantic`) /
  sd(cell_ua$Sycophantic - cell_ua$`Anti-sycophantic`)
c(M_syco = mean(cell_ua$Sycophantic), M_anti = mean(cell_ua$`Anti-sycophantic`))

## Obsequiousness
cell_ob <- dat %>%
  group_by(feedback_id, condition) %>%
  summarise(m = mean(Obsequiousness, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = m) %>%
  arrange(feedback_id)
t.test(cell_ob$Sycophantic, cell_ob$`Anti-sycophantic`, paired = TRUE)
mean(cell_ob$Sycophantic - cell_ob$`Anti-sycophantic`) /
  sd(cell_ob$Sycophantic - cell_ob$`Anti-sycophantic`)
c(M_syco = mean(cell_ob$Sycophantic), M_anti = mean(cell_ob$`Anti-sycophantic`))

## Excitement
cell_ex <- dat %>%
  group_by(feedback_id, condition) %>%
  summarise(m = mean(Excitement, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = m) %>%
  arrange(feedback_id)
t.test(cell_ex$Sycophantic, cell_ex$`Anti-sycophantic`, paired = TRUE)
mean(cell_ex$Sycophantic - cell_ex$`Anti-sycophantic`) /
  sd(cell_ex$Sycophantic - cell_ex$`Anti-sycophantic`)
c(M_syco = mean(cell_ex$Sycophantic), M_anti = mean(cell_ex$`Anti-sycophantic`))

# ============================================================================
# SSS Validation study and PCA Exp 2 correlation 
# (recover which feedbacks each exp 2 participant saw, average the SSS of those feedbacks, run the PCA, correlate)
# ============================================================================
COL_SYC <- "#A8E6A3"   # sycophantic 
COL_CON <- "#F4A6A6"   # anti-sycophantic 

## SSS per feedback x condition (mean over raters)
human_cell <- dat %>%
  group_by(feedback_id, condition) %>%
  summarise(H_Overall_sycophancy = mean(Overall_sycophancy), .groups = "drop") %>%
  rename(fid_hll = feedback_id, cond = condition) %>%
  mutate(cond = as.character(cond))

## Exp 2 causal participant data + feedback coding
cz        <- read_csv("data/data_experiment_2_caus_vd - dat.csv", show_col_types = FALSE)
code_caus <- read_csv("data/Causal_Feedback_Coding_Exp2_vd - dat.csv", show_col_types = FALSE)

cz <- cz %>%
  mutate(pid = row_number(),
         cond = case_when(condition == "sycophantic"     ~ "Sycophantic",
                          condition == "non_sycophantic" ~ "Anti-sycophantic",
                          TRUE ~ NA_character_))

## Lookups: (question, answer value) -> causal feedback_id -> human id
code_caus  <- code_caus %>% mutate(q = as.integer(str_extract(as.character(`Q number`), "\\d+")))
val_lookup <- code_caus %>% transmute(q, val = as.integer(Coding), fid_orig = feedback_id)
crosswalk  <- code_caus %>% transmute(fid_orig = feedback_id,
                                      fid_hll  = coalesce(`In LLM_HumanRaters`, feedback_id))

## Recover which feedback each participant saw, attach SSS, average per participant
long <- cz %>%
  select(pid, cond, all_of(paste0("Q", 1:10, "_pre"))) %>%
  pivot_longer(starts_with("Q"), names_to = "q", values_to = "val") %>%
  mutate(q = as.integer(str_extract(q, "\\d+")), val = as.integer(val)) %>%
  left_join(val_lookup, by = c("q", "val")) %>%
  left_join(crosswalk,  by = "fid_orig") %>%
  left_join(human_cell, by = c("fid_hll", "cond"))

part_scores <- long %>%
  group_by(pid) %>%
  summarise(part_H_Overall = mean(H_Overall_sycophancy, na.rm = TRUE), .groups = "drop")
cz <- cz %>% left_join(part_scores, by = "pid")

## PCA on the five tone ratings -> per-participant sycophancy_score
items5 <- c("empathetic_1", "polite_1", "helpful_1", "blunt_1", "confrontational_1")
tone5  <- cz %>%
  select(all_of(items5)) %>%
  na.omit() %>%
  mutate(across(everything(), as.numeric))

# Sampling adequacy
KMO(tone5)
cortest.bartlett(cor(tone5), n = nrow(tone5))

# Parallel analysis (simulated + permuted nulls, 95th percentile)
n <- nrow(tone5); p <- ncol(tone5)
obs_eigen <- eigen(cor(tone5))$values
set.seed(2026)
n_iter  <- 1000
eig_sim <- matrix(NA_real_, n_iter, p)
eig_res <- matrix(NA_real_, n_iter, p)
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

# Loadings and variance explained
pca_result <- principal(tone5, nfactors = 1, rotate = "none")
print(pca_result$loadings, cutoff = 0, digits = 3)
pca_result$Vaccounted

# Component scores 
pc_fit     <- prcomp(tone5, center = TRUE, scale. = TRUE)
pc1_scores <- pc_fit$x[, 1]
if (cor(pc1_scores, tone5$empathetic_1) < 0) pc1_scores <- -pc1_scores
cz$sycophancy_score <- NA_real_
cc <- complete.cases(cz[, items5])
cz$sycophancy_score[cc] <- pc1_scores

## Pearson correlation 
cor.test(cz$sycophancy_score, cz$part_H_Overall)

# ============================================================================
#Figures
# ============================================================================

## ---- Fig 4a: PCA loadings (PC1) ----
tone_labels_map <- c(empathetic_1 = "Empathy", polite_1 = "Politeness",
                     helpful_1 = "Helpfulness", blunt_1 = "Bluntness",
                     confrontational_1 = "Confrontation")
L <- as.numeric(unclass(pca_result$loadings)[, 1])
names(L) <- rownames(unclass(pca_result$loadings))
if (L["empathetic_1"] < 0) L <- -L
loadings_df <- data.frame(Variable = names(L), Loading = L) %>%
  mutate(Variable = recode(Variable, !!!tone_labels_map),
         Variable = factor(Variable, levels = c("Empathy","Politeness","Helpfulness",
                                                "Bluntness","Confrontation")))

ggplot(loadings_df, aes(x = Variable, y = Loading)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = sprintf("%.2f", Loading),
                vjust = ifelse(Loading >= 0, -0.5, 1.3)), size = 4) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = NULL, y = "Factor Loadings") +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(size = 28, angle = 35, hjust = 1, face = "bold"),
        axis.title.y = element_text(size = 28, face = "bold"),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

## ---- Fig 4b: SSS violin (total + 3 subscales, dodged by condition) ----
# Each dot = one rating (per the caption). Uses per-rating scores in `dat`.
## ---- Fig 4b: SSS violin (total + 3 subscales, dodged by condition) ----
# Each dot = one rating (per the caption). Uses per-rating scores in `dat`.
sss_plot_df <- dat %>%
  transmute(condition,
            Sycophancy             = Overall_sycophancy,
            Obsequiousness         = Obsequiousness,
            Excitement             = Excitement,
            `Uncritical Agreement` = Uncritical_Agreement) %>%
  pivot_longer(-condition, names_to = "scale", values_to = "rating") %>%
  mutate(scale = factor(scale, levels = c("Sycophancy","Obsequiousness",
                                          "Excitement","Uncritical Agreement"),
                        labels = c("Social\nSycophancy Scale","Obsequiousness",
                                   "Excitement","Uncritical\nAgreement")),
         condition = factor(ifelse(condition == "Sycophantic",
                                   "Sycophantic AI tone","Anti-sycophantic AI tone"),
                            levels = c("Sycophantic AI tone","Anti-sycophantic AI tone")))
ggplot(sss_plot_df, aes(x = scale, y = rating, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = "grey40", linewidth = 0.3,
              position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.9, seed = 2026),
              alpha = 0.35, size = 0.9, colour = "grey20") +
  stat_summary(aes(group = condition), fun.data = mean_se, geom = "errorbar",
               width = 0.12, linewidth = 0.6, position = position_dodge(width = 0.9)) +
  stat_summary(aes(group = condition), fun = mean, geom = "point", shape = 23,
               size = 3, fill = "white", colour = "black", stroke = 0.8,
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Sycophantic AI tone" = COL_SYC,
                               "Anti-sycophantic AI tone" = COL_CON)) +
  scale_y_continuous(breaks = 1:5) +
  coord_cartesian(ylim = c(0.8, 5.2)) +
  labs(x = NULL, y = "Ratings", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top", axis.text.x = element_text(face = "bold"))


## ---- Fig 4c: PCA score by condition (violin) ----
pca_plot_df <- cz %>%
  filter(!is.na(sycophancy_score), !is.na(cond)) %>%
  mutate(Condition = factor(
    recode(cond, "Sycophantic" = "Sycophantic\nAI tone",
           "Anti-sycophantic" = "Anti-sycophantic\nAI tone"),
    levels = c("Sycophantic\nAI tone", "Anti-sycophantic\nAI tone")))

ggplot(pca_plot_df, aes(x = Condition, y = sycophancy_score, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_jitter(height = 0, width = 0.08, alpha = 0.35, size = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, fill = "white") +
  scale_fill_manual(values = c("Sycophantic\nAI tone" = COL_SYC,
                               "Anti-sycophantic\nAI tone" = COL_CON)) + 
  labs(x = NULL, y = "PCA score") + 
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"))

