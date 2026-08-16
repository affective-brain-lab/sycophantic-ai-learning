# AI’s sycophantic tone makes users feel smarter without becoming smarter

Data and code from manuscript: Duettmann, A. & Sharot, T. (in prep.). AI’s sycophantic tone makes users feel smarter without becoming smarter

## Overview
Two experiments (N = 848) examined how AI’s sycophantic tone influences subjective and objective aspects of human learning. Experiment 1 (N = 288) used a correlational design in which participants interacted with AI responses with natural variation in perceived sycophantic tone measured via participants’ ratings and related to learning outcomes. Experiment 2 (N = 256) used a between-subjects design in which participants interacted either with sycophantic or anti-sycophantic AI responses, testing the relationships causally. Both experiments measured subjective aspects of learning (mood, perceived competence, confidence and the desire to learn) alongside an objective transfer test of knowledge acquisition. Control groups for each experiment completed the studies without AI interaction. A separate validation study with human raters confirmed that the AI responses produced perceived social sycophancy using an established social sycophancy scale.

## Structure

- `data/` — anonymised participant data (CSV)
  - `data_experiment_1_corr_Main_vf - dat.csv` 
  - `data_experiment_1_corr_wControl_vf - dat.csv`
  - `data_experiment_2_caus_vd - dat.csv`
  - `data_experiment_2_caus_wControl_vd - dat.csv` 
  - `Causal_Feedback_Coding_Exp2_vd - dat.csv` — Experiment II feedback coding
  - `data_vadlidation_study_vd - dat.csv` — validation study ratings
- `scripts/`
  - `Exp1_Corr_Script_vd.R` — Experiment I (correlational)
  - `Exp2_Caus_Script_vd.R` — Experiment II (causal)
  - `Validation_Study_Script_vd.R` — validation study 
- `README.md`


## Requirements

R (≥ 4.x) with: `tidyverse`, `psych`, `BayesFactor`, `mediation`, `lmtest`,
`sandwich`, `car`, `gt`, `ggtext`, `conflicted`.
