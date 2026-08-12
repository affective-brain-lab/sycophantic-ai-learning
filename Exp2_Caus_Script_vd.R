# ============================================================================
# EXPERIMENT II (causal)
#Excl OLS assumptions + HC3 supp material table
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
conflict_prefer("summarise","dplyr")
conflict_prefer("mutate",  "dplyr")
conflict_prefer("mediate", "mediation")

# ============================================================================
# 1. Data prep 
# ============================================================================
df <- read_csv("/Users/a.duett/Downloads/Causal_Phil_03_positive_neutral_May+14,+2026_14.20 - DAT_all_pos_neg_comb_final - DAT_all_pos_neg_comb_f_ACs.csv")

pre_cols        <- paste0("Q", 1:10, "_pre")            # initial exam raw responses
post_cols       <- paste0("Q", 1:10, "_post")           # final exam raw responses
confi_pre_cols  <- paste0("Q", 1:10, "_confi_pre_1")    # per-question confidence, initial
confi_post_cols <- paste0("Q", 1:10, "_confi_1")        # per-question confidence, final
learn_cols      <- paste0("Q", 1:10, "_learn_1")        # per-question desire to learn

df <- df %>%
  mutate(across(all_of(c(pre_cols, post_cols,
                         confi_pre_cols, confi_post_cols, learn_cols)),
                ~ as.numeric(.)))

df <- df %>%
  mutate(
    # Initial exam accuracy: Q1/Q2 correct == 3; Q3-Q10 correct == 1
    correct_Q1_pre  = as.integer(Q1_pre  == 3),
    correct_Q2_pre  = as.integer(Q2_pre  == 3),
    correct_Q3_pre  = as.integer(Q3_pre  == 1),
    correct_Q4_pre  = as.integer(Q4_pre  == 1),
    correct_Q5_pre  = as.integer(Q5_pre  == 1),
    correct_Q6_pre  = as.integer(Q6_pre  == 1),
    correct_Q7_pre  = as.integer(Q7_pre  == 1),
    correct_Q8_pre  = as.integer(Q8_pre  == 1),
    correct_Q9_pre  = as.integer(Q9_pre  == 1),
    correct_Q10_pre = as.integer(Q10_pre == 1),
    accuracy_pre  = rowSums(across(starts_with("correct_Q") & ends_with("_pre")),
                            na.rm = TRUE) / 10,
    # Final exam accuracy: correct == 1 for all 10
    accuracy_post = rowSums(across(all_of(post_cols)) == 1, na.rm = TRUE) / 10,
    avg_confi_pre  = rowMeans(across(all_of(confi_pre_cols)),  na.rm = TRUE),
    avg_confi_post = rowMeans(across(all_of(confi_post_cols)), na.rm = TRUE),
    avg_learn      = rowMeans(across(all_of(learn_cols)),      na.rm = TRUE)
  ) %>%
  select(-starts_with("correct_Q"))

# ============================================================================
# 2. Variables
# ============================================================================
df <- df %>%
  rename(
    happy_post            = happy_post_1,
    happy_pre             = happy_pre_1,
    perceived_comp_post   = understanding_post_1,
    perceived_comp_pre    = understanding_pre_1,
    desire_learn_phil     = avg_learn,          # studied domain (philosophy)
    desire_learn_art      = interest_post_1,    # novel domain (art history)
    est_mark_initialexam  = confidence_1_1,     # % correct estimate, initial exam
    est_mark_transferexam = Confidence_2_1      # % correct estimate, final exam
  ) %>%
  mutate(
    condition = case_when(
      condition == "syco"     ~ "sycophantic",
      condition == "non-syco" ~ "anti_sycophantic",
      TRUE                    ~ NA_character_
    ),
    condition = factor(condition, levels = c("anti_sycophantic", "sycophantic"))
  )

df %>% count(condition)   # sanity check

# ============================================================================
# 3. PCA on the five AI-tone ratings 
# ============================================================================
items5 <- c("empathetic_1", "polite_1", "helpful_1", "blunt_1", "confrontational_1")
tone5  <- df %>%
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

#Loadings and variance explained
pca_result <- principal(tone5, nfactors = 1, rotate = "none")
print(pca_result$loadings, cutoff = 0, digits = 3)
pca_result$Vaccounted

# Component scores 
pc_fit     <- prcomp(tone5, center = TRUE, scale. = TRUE)
pc1_scores <- pc_fit$x[, 1]
if (cor(pc1_scores, tone5$empathetic_1) < 0) pc1_scores <- -pc1_scores
df$sycophancy_score <- NA_real_
cc <- complete.cases(df[, items5])
df$sycophancy_score[cc] <- pc1_scores

# Does the PCA score differ by condition?
df %>% group_by(condition) %>%
  summarise(n = n(), M = mean(sycophancy_score, na.rm = TRUE),
            SD = sd(sycophancy_score, na.rm = TRUE))
t.test(sycophancy_score ~ condition, data = df, var.equal = FALSE)   

# ============================================================================
# 4. Z-score continuous outcomes / covariates 
# ============================================================================
z_vars <- c("accuracy_pre", "accuracy_post",
            "happy_pre", "happy_post",
            "perceived_comp_pre", "perceived_comp_post",
            "avg_confi_pre", "avg_confi_post",
            "desire_learn_phil", "desire_learn_art",
            "est_mark_initialexam", "est_mark_transferexam")

df_z <- df %>%
  mutate(across(all_of(z_vars), ~ as.numeric(scale(.))))

# ============================================================================
# 5.MODELS
# ============================================================================

## S3.1 Mood
model_S3_1a <- lm(happy_post ~ condition + happy_pre, data = df_z)
summary(model_S3_1a)
model_S3_1b <- lm(happy_post ~ condition + happy_pre + accuracy_pre + accuracy_post, data = df_z)
summary(model_S3_1b)

## S3.2 Perceived competence
model_S3_2a <- lm(perceived_comp_post ~ condition + perceived_comp_pre, data = df_z)
summary(model_S3_2a)
model_S3_2b <- lm(perceived_comp_post ~ condition + perceived_comp_pre + accuracy_pre + accuracy_post, data = df_z)
summary(model_S3_2b)

## S3.3 Desire to learn (philosophy)
model_S3_3a <- lm(desire_learn_phil ~ condition, data = df_z)
summary(model_S3_3a)
model_S3_3b <- lm(desire_learn_phil ~ condition + accuracy_pre + accuracy_post, data = df_z)
summary(model_S3_3b)

## S3.4 Desire to learn (art history / novel domain)
model_S3_4a <- lm(desire_learn_art ~ condition, data = df_z)
summary(model_S3_4a)
model_S3_4b <- lm(desire_learn_art ~ condition + accuracy_pre + accuracy_post, data = df_z)
summary(model_S3_4b)

## S3.5 Confidence
model_S3_5a <- lm(avg_confi_post ~ condition + avg_confi_pre, data = df_z)
summary(model_S3_5a)
model_S3_5b <- lm(avg_confi_post ~ condition + avg_confi_pre + accuracy_pre + accuracy_post, data = df_z)
summary(model_S3_5b)

## S3.6 Estimated final mark
model_S3_6a <- lm(est_mark_transferexam ~ condition + est_mark_initialexam, data = df_z)
summary(model_S3_6a)
model_S3_6b <- lm(est_mark_transferexam ~ condition + est_mark_initialexam + accuracy_pre + accuracy_post, data = df_z)
summary(model_S3_6b)

## S3.7 Actual competence 
model_S3_7a <- lm(accuracy_post ~ condition, data = df_z)
summary(model_S3_7a)
model_S3_7b <- lm(accuracy_post ~ condition + accuracy_pre, data = df_z)
summary(model_S3_7b)

## S3.8 Final exam performance: AI conditions vs. no-feedback control

## Dataset with control group:
df_ctrl_raw <- read_csv("/Users/a.duett/Downloads/data_experiment_2_causal_wControl_vd - dat.csv")

df_ctrl_z <- df_ctrl_raw %>%
  mutate(
    condition = case_when(
      condition == "sycophantic"     ~ "sycophantic",
      condition == "non_sycophantic" ~ "anti_sycophantic",
      condition == "control"         ~ "control",
      TRUE                           ~ NA_character_
    ),
    condition = factor(condition, levels = c("control", "sycophantic", "anti_sycophantic"))
  ) %>%
  mutate(across(c(accuracy_pre, accuracy_post), ~ as.numeric(scale(.))))   # z across all 3

## --- Sycophantic vs. no-feedback control ---
# control = reference
df_syc_base <- df_ctrl_z %>%
  filter(condition %in% c("control", "sycophantic")) %>%
  mutate(condition = factor(condition, levels = c("control", "sycophantic")))

model_S3_8a_syc <- lm(accuracy_post ~ condition, data = df_syc_base)
summary(model_S3_8a_syc)
model_S3_8b_syc <- lm(accuracy_post ~ condition + accuracy_pre, data = df_syc_base)
summary(model_S3_8b_syc)

## --- Anti-sycophantic vs. no-feedback control ---
# control = reference
df_anti_base <- df_ctrl_z %>%
  filter(condition %in% c("control", "anti_sycophantic")) %>%
  mutate(condition = factor(condition, levels = c("control", "anti_sycophantic")))

model_S3_8a_anti <- lm(accuracy_post ~ condition, data = df_anti_base)
summary(model_S3_8a_anti)
model_S3_8b_anti <- lm(accuracy_post ~ condition + accuracy_pre, data = df_anti_base)
summary(model_S3_8b_anti)

## initial exam performance 
model_S3_8_pre_syc  <- lm(accuracy_pre ~ condition, data = df_syc_base)
summary(model_S3_8_pre_syc)
model_S3_8_pre_anti <- lm(accuracy_pre ~ condition, data = df_anti_base)
summary(model_S3_8_pre_anti)

# ============================================================================
# 6. BAYES FACTORS  (BF01 > 1)
# ============================================================================
df_bf <- df_z %>% as.data.frame() %>%
  mutate(condition = droplevels(factor(condition)))

## Confidence 
df_bf_c <- df_bf[complete.cases(df_bf[, c("avg_confi_post","condition","avg_confi_pre")]), ]
set.seed(2026)
bf_full_confi  <- lmBF(avg_confi_post ~ condition + avg_confi_pre, data = df_bf_c)
bf_reduc_confi <- lmBF(avg_confi_post ~ avg_confi_pre,             data = df_bf_c)
bf01_confi     <- bf_reduc_confi / bf_full_confi
bf01_confi

## Estimated final mark 
df_bf_e <- df_bf[complete.cases(df_bf[, c("est_mark_transferexam","condition","est_mark_initialexam")]), ]
set.seed(2026)
bf_full_est  <- lmBF(est_mark_transferexam ~ condition + est_mark_initialexam, data = df_bf_e)
bf_reduc_est <- lmBF(est_mark_transferexam ~ est_mark_initialexam,             data = df_bf_e)
bf01_est     <- bf_reduc_est / bf_full_est
bf01_est

## Actual competence 
df_bf_a <- df_bf[complete.cases(df_bf[, c("accuracy_post","condition","accuracy_pre")]), ]
set.seed(2026)
bf_acc_nc  <- anovaBF(accuracy_post ~ condition, data = df_bf_a)
bf01_acc_nc <- 1 / bf_acc_nc
bf01_acc_nc

set.seed(2026)
bf_full_acc  <- lmBF(accuracy_post ~ condition + accuracy_pre, data = df_bf_a)
bf_reduc_acc <- lmBF(accuracy_post ~ accuracy_pre,             data = df_bf_a)
bf01_acc_ctrl <- bf_reduc_acc / bf_full_acc
bf01_acc_ctrl

# ============================================================================
# 7. MEDIATION MODELS 
# ============================================================================

## S4.1 Perceived competence
# a) 
med_m_S4_1a <- lm(happy_post ~ condition + happy_pre, data = df_z)
out_m_S4_1a <- lm(perceived_comp_post ~ condition + happy_post + happy_pre +
                    perceived_comp_pre, data = df_z)
set.seed(2026)
med_S4_1a <- mediate(med_m_S4_1a, out_m_S4_1a,
                     treat = "condition", mediator = "happy_post",
                     control.value = "anti_sycophantic", treat.value = "sycophantic",
                     boot = TRUE, sims = 1000)
summary(med_S4_1a)

# b) 
med_m_S4_1b <- lm(happy_post ~ condition + happy_pre + accuracy_pre + accuracy_post, data = df_z)
out_m_S4_1b <- lm(perceived_comp_post ~ condition + happy_post + happy_pre +
                    perceived_comp_pre + accuracy_pre + accuracy_post, data = df_z)
set.seed(2026)
med_S4_1b <- mediate(med_m_S4_1b, out_m_S4_1b,
                     treat = "condition", mediator = "happy_post",
                     control.value = "anti_sycophantic", treat.value = "sycophantic",
                     boot = TRUE, sims = 1000)
summary(med_S4_1b)

## S4.2 Desire to learn (philosophy)
# a) 
med_m_S4_2a <- lm(happy_post ~ condition + happy_pre, data = df_z)
out_m_S4_2a <- lm(desire_learn_phil ~ condition + happy_post + happy_pre, data = df_z)
set.seed(2026)
med_S4_2a <- mediate(med_m_S4_2a, out_m_S4_2a,
                     treat = "condition", mediator = "happy_post",
                     control.value = "anti_sycophantic", treat.value = "sycophantic",
                     boot = TRUE, sims = 1000)
summary(med_S4_2a)

# b) 
med_m_S4_2b <- lm(happy_post ~ condition + happy_pre + accuracy_pre + accuracy_post, data = df_z)
out_m_S4_2b <- lm(desire_learn_phil ~ condition + happy_post + happy_pre +
                    accuracy_pre + accuracy_post, data = df_z)
set.seed(2026)
med_S4_2b <- mediate(med_m_S4_2b, out_m_S4_2b,
                     treat = "condition", mediator = "happy_post",
                     control.value = "anti_sycophantic", treat.value = "sycophantic",
                     boot = TRUE, sims = 1000)
summary(med_S4_2b)

## S4.3 Desire to learn (art history)
# a) 
med_m_S4_3a <- lm(happy_post ~ condition + happy_pre, data = df_z)
out_m_S4_3a <- lm(desire_learn_art ~ condition + happy_post + happy_pre, data = df_z)
set.seed(2026)
med_S4_3a <- mediate(med_m_S4_3a, out_m_S4_3a,
                     treat = "condition", mediator = "happy_post",
                     control.value = "anti_sycophantic", treat.value = "sycophantic",
                     boot = TRUE, sims = 1000)
summary(med_S4_3a)

# b) 
med_m_S4_3b <- lm(happy_post ~ condition + happy_pre + accuracy_pre + accuracy_post, data = df_z)
out_m_S4_3b <- lm(desire_learn_art ~ condition + happy_post + happy_pre +
                    accuracy_pre + accuracy_post, data = df_z)
set.seed(2026)
med_S4_3b <- mediate(med_m_S4_3b, out_m_S4_3b,
                     treat = "condition", mediator = "happy_post",
                     control.value = "anti_sycophantic", treat.value = "sycophantic",
                     boot = TRUE, sims = 1000)
summary(med_S4_3b)



# ============================================================================
# AGE robustness check 
# ============================================================================
age_col <- "Age"                      
summary(as.numeric(df$Age))

z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
df$age   <- suppressWarnings(as.numeric(df[[age_col]]))
df_z$age <- z(df$age)                 # df_z is df + mutate only -> rows align

star <- function(p) ifelse(is.na(p), "",
                           ifelse(p < .001, "***", ifelse(p < .01, "**",
                                                          ifelse(p < .05, "*", ""))))

# fit base and base+age on the SAME complete-case subset
focal_row <- function(label, outcome, rhs, data, focal) {
  vars <- unique(c(outcome, all.vars(as.formula(paste("~", rhs))), "age"))
  d  <- data[stats::complete.cases(data[, vars]), ]
  m0 <- lm(as.formula(paste(outcome, "~", rhs)), data = d)
  m1 <- lm(as.formula(paste(outcome, "~", rhs, "+ age")), data = d)
  g  <- function(m, term) {
    s <- summary(m)$coefficients
    if (!term %in% rownames(s)) return(c(b = NA, p = NA))
    t <- s[term, "t value"]
    c(b = s[term, "Estimate"], p = s[term, "Pr(>|t|)"])
  }
  a0 <- g(m0, focal); a1 <- g(m1, focal)
  sm <- summary(m1)$coefficients
  age_p <- if ("age" %in% rownames(sm)) sm["age", "Pr(>|t|)"] else NA
  data.frame(
    model    = label, N = nrow(d),
    b_noAge  = round(a0["b"], 3), p_noAge = round(a0["p"], 4), s_noAge = star(a0["p"]),
    b_age    = round(a1["b"], 3), p_age   = round(a1["p"], 4), s_age   = star(a1["p"]),
    dp       = round(a1["p"] - a0["p"], 4),
    flipped  = (a0["p"] < .05) != (a1["p"] < .05),
    age_p    = round(age_p, 4),
    row.names = NULL)
}

# ---- 1. MAIN MODELS (df_z: sycophantic vs anti-sycophantic) -----------------
S <- "conditionsycophantic"
main_compare <- rbind(
  focal_row("S3.1a Mood",            "happy_post",            "condition + happy_pre", df_z, S),
  focal_row("S3.1b Mood",            "happy_post",            "condition + happy_pre + accuracy_pre + accuracy_post", df_z, S),
  focal_row("S3.2a Perceived comp",  "perceived_comp_post",   "condition + perceived_comp_pre", df_z, S),
  focal_row("S3.2b Perceived comp",  "perceived_comp_post",   "condition + perceived_comp_pre + accuracy_pre + accuracy_post", df_z, S),
  focal_row("S3.3a Desire (phil)",   "desire_learn_phil",     "condition", df_z, S),
  focal_row("S3.3b Desire (phil)",   "desire_learn_phil",     "condition + accuracy_pre + accuracy_post", df_z, S),
  focal_row("S3.4a Desire (art)",    "desire_learn_art",      "condition", df_z, S),
  focal_row("S3.4b Desire (art)",    "desire_learn_art",      "condition + accuracy_pre + accuracy_post", df_z, S),
  focal_row("S3.5a Confidence",      "avg_confi_post",        "condition + avg_confi_pre", df_z, S),
  focal_row("S3.5b Confidence",      "avg_confi_post",        "condition + avg_confi_pre + accuracy_pre + accuracy_post", df_z, S),
  focal_row("S3.6a Est. mark",       "est_mark_transferexam", "condition + est_mark_initialexam", df_z, S),
  focal_row("S3.6b Est. mark",       "est_mark_transferexam", "condition + est_mark_initialexam + accuracy_pre + accuracy_post", df_z, S),
  focal_row("S3.7a Actual comp",     "accuracy_post",         "condition", df_z, S),
  focal_row("S3.7b Actual comp",     "accuracy_post",         "condition + accuracy_pre", df_z, S)
)

# ---- 1b. S3.8 control comparisons (only if age is in the control file) ------
if (exists("df_ctrl_raw") && age_col %in% names(df_ctrl_raw)) {
  df_ctrl_z$age <- z(suppressWarnings(as.numeric(df_ctrl_raw[[age_col]])))
  df_syc_base  <- df_ctrl_z %>% filter(condition %in% c("control","sycophantic")) %>%
    mutate(condition = factor(condition, levels = c("control","sycophantic")))
  df_anti_base <- df_ctrl_z %>% filter(condition %in% c("control","anti_sycophantic")) %>%
    mutate(condition = factor(condition, levels = c("control","anti_sycophantic")))
  main_compare <- rbind(main_compare,
                        focal_row("S3.8a Syc vs ctrl",  "accuracy_post", "condition",                df_syc_base,  "conditionsycophantic"),
                        focal_row("S3.8b Syc vs ctrl",  "accuracy_post", "condition + accuracy_pre", df_syc_base,  "conditionsycophantic"),
                        focal_row("S3.8a Anti vs ctrl", "accuracy_post", "condition",                df_anti_base, "conditionanti_sycophantic"),
                        focal_row("S3.8b Anti vs ctrl", "accuracy_post", "condition + accuracy_pre", df_anti_base, "conditionanti_sycophantic"))
} else message("Age not found in control-group file - S3.8 age check skipped.")

print(main_compare, row.names = FALSE)
# write_csv(main_compare, "age_check_main_models.csv")

# ---- 2. MEDIATION MODELS ----------------------------------------------------
med_row <- function(m) data.frame(
  ACME = m$d.avg, ACME_p = m$d.avg.p,
  ADE  = m$z.avg, ADE_p  = m$z.avg.p,
  prop = m$n.avg, prop_p = m$n.avg.p)


run_med <- function(med_rhs, out_rhs, outcome, data, add_age = FALSE) {
  if (add_age) { med_rhs <- paste(med_rhs, "+ age"); out_rhs <- paste(out_rhs, "+ age") }
  med_form <- as.formula(paste("happy_post ~", med_rhs))
  out_form <- as.formula(paste(outcome, "~", out_rhs))
  # force the formula to be evaluated and stored in the call, not a variable ref
  mm <- do.call("lm", list(formula = med_form, data = quote(data)))
  om <- do.call("lm", list(formula = out_form, data = quote(data)))
  environment(mm$terms) <- environment()
  environment(om$terms) <- environment()
  set.seed(2026)
  mediate(mm, om, treat = "condition", mediator = "happy_post",
          control.value = "anti_sycophantic", treat.value = "sycophantic",
          boot = TRUE, sims = 1000)
}

compare_med <- function(label, outcome, med_rhs, out_rhs, data = df_z) {
  vars <- unique(c("happy_post", outcome, "condition", "age",
                   all.vars(as.formula(paste("~", med_rhs))),
                   all.vars(as.formula(paste("~", out_rhs)))))
  d  <- data[stats::complete.cases(data[, vars]), ]
  m0 <- run_med(med_rhs, out_rhs, outcome, d, FALSE)
  m1 <- run_med(med_rhs, out_rhs, outcome, d, TRUE)
  r  <- rbind(cbind(model = paste(label, "| no age"), med_row(m0)),
              cbind(model = paste(label, "| + age"),  med_row(m1)))
  r[, -1] <- round(r[, -1], 4); r
}

med_compare <- rbind(
  compare_med("S4.1a Perceived comp", "perceived_comp_post",
              "condition + happy_pre",
              "condition + happy_post + happy_pre + perceived_comp_pre"),
  compare_med("S4.1b Perceived comp", "perceived_comp_post",
              "condition + happy_pre + accuracy_pre + accuracy_post",
              "condition + happy_post + happy_pre + perceived_comp_pre + accuracy_pre + accuracy_post"),
  compare_med("S4.2a Desire (phil)",  "desire_learn_phil",
              "condition + happy_pre",
              "condition + happy_post + happy_pre"),
  compare_med("S4.2b Desire (phil)",  "desire_learn_phil",
              "condition + happy_pre + accuracy_pre + accuracy_post",
              "condition + happy_post + happy_pre + accuracy_pre + accuracy_post"),
  compare_med("S4.3a Desire (art)",   "desire_learn_art",
              "condition + happy_pre",
              "condition + happy_post + happy_pre"),
  compare_med("S4.3b Desire (art)",   "desire_learn_art",
              "condition + happy_pre + accuracy_pre + accuracy_post",
              "condition + happy_post + happy_pre + accuracy_pre + accuracy_post")
)

print(med_compare, row.names = FALSE)

# ============================================================================
# 10. FIGURES 
# ============================================================================
col_syco <- "#A8E6A3"   # sycophantic 
col_anti <- "#F4A6A6"   # anti-sycophantic 
col_base <- "#CCCCCC"   # no-feedback control

## ---- Fig 4a: PCA loadings (PC1) ----
tone_labels_map <- c(empathetic_1 = "Empathy", polite_1 = "Politeness",
                     helpful_1 = "Helpfulness", blunt_1 = "Bluntness",
                     confrontational_1 = "Confrontation")
tone_levels <- c("Empathy", "Politeness", "Helpfulness", "Bluntness", "Confrontation")

L <- as.numeric(unclass(pca_result$loadings)[, 1])
names(L) <- rownames(unclass(pca_result$loadings))
if (L["empathetic_1"] < 0) L <- -L
loadings_df <- data.frame(Variable = names(L), Loading = L) %>%
  mutate(Variable = recode(Variable, !!!tone_labels_map),
         Variable = factor(Variable, levels = tone_levels))

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

## ---- Fig 4c: PCA score by condition (violin) ----
df_pca_plot <- df %>%
  filter(!is.na(sycophancy_score)) %>%
  mutate(Condition = factor(
    recode(as.character(condition),
           "sycophantic" = "Sycophantic\nAI tone",
           "anti_sycophantic" = "Anti-sycophantic\nAI tone"),
    levels = c("Sycophantic\nAI tone", "Anti-sycophantic\nAI tone")))

ggplot(df_pca_plot, aes(x = Condition, y = sycophancy_score, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_jitter(height = 0, width = 0.08, alpha = 0.35, size = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, fill = "white") +
  scale_fill_manual(values = c("Sycophantic\nAI tone" = col_syco,
                               "Anti-sycophantic\nAI tone" = col_anti)) +
  labs(x = NULL, y = "Perceived Sycophantic Tone (PCA score)") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"))

## ---- Fig 5: subjective outcomes (sycophantic vs. anti-sycophantic) ----
df_plot2 <- df_z %>%
  mutate(Condition = factor(
    recode(as.character(condition),
           "sycophantic" = "Sycophantic tone",
           "anti_sycophantic" = "Anti-sycophantic tone"),
    levels = c("Sycophantic tone", "Anti-sycophantic tone")))

cond_colors2 <- c("Sycophantic tone" = col_syco, "Anti-sycophantic tone" = col_anti)

theme_bar <- theme_classic(base_size = 14) +
  theme(axis.title.y = element_markdown(size = 26, margin = margin(r = 4)),
        axis.text.x  = element_text(size = 18, face = "bold"),
        legend.position = "none")

# Fig 5a - Mood
ggplot(df_plot2, aes(x = Condition, y = happy_post, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.85, width = 0.65) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = cond_colors2) +
  labs(x = NULL, y = "**Mood**<br><span style='font-size:16pt'>(post AI-interaction, z-scored)</span>") +
  theme_bar

# Fig 5b - Perceived competence
ggplot(df_plot2, aes(x = Condition, y = perceived_comp_post, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.85, width = 0.65) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = cond_colors2) +
  labs(x = NULL, y = "**Perceived competence**<br><span style='font-size:16pt'>(post AI-interaction, z-scored)</span>") +
  theme_bar

# Fig 5c - Desire to learn (philosophy + art history, grouped)
df_desire_long <- df_plot2 %>%
  select(Condition, desire_learn_phil, desire_learn_art) %>%
  pivot_longer(c(desire_learn_phil, desire_learn_art),
               names_to = "Domain", values_to = "desire_z") %>%
  mutate(Domain = factor(recode(Domain,
                                "desire_learn_phil" = "Studied domain\n(Philosophy)",
                                "desire_learn_art"  = "Novel domain\n(Art History)"),
                         levels = c("Studied domain\n(Philosophy)", "Novel domain\n(Art History)")))

ggplot(df_desire_long, aes(x = Domain, y = desire_z, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.72),
               color = "black", alpha = 0.85, width = 0.65) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.72),
               width = 0.18, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = cond_colors2, name = NULL) +
  scale_x_discrete(labels = function(x) gsub("\n", "<br>", x)) +
  labs(x = NULL, y = "**Desire to learn**<br><span style='font-size:16pt'>(post AI-interaction, z-scored)</span>") +
  theme_bar %+replace%
  theme(axis.text.x = element_markdown(size = 13, face = "bold"),
        legend.position = "right", legend.text = element_text(size = 11))

## ---- Fig 6a: actual competence ----
df_plot3 <- df_ctrl_z %>%
  mutate(Condition = factor(
    recode(as.character(condition),
           "sycophantic" = "Sycophantic tone",
           "anti_sycophantic" = "Anti-sycophantic tone",
           "baseline" = "Without AI feedback"),
    levels = c("Sycophantic tone", "Anti-sycophantic tone", "Without AI feedback")))

cond_colors3 <- c("Sycophantic tone" = col_syco, "Anti-sycophantic tone" = col_anti,
                  "Without AI feedback" = col_base)

ggplot(df_plot3, aes(x = Condition, y = accuracy_post, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.85, width = 0.65) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = cond_colors3) +
  labs(x = NULL, y = "**Competence**<br><span style='font-size:16pt'>(final exam score, z-scored)</span>") +
  theme_bar %+replace%
  theme(axis.text.x = element_text(size = 15, face = "bold", angle = 20, hjust = 1))



