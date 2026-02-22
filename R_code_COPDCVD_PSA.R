# COPD-CVD Markov Model: Probabilistic Sensitivity Analysis
# 19-state model with 1,000 PSA iterations

library(tidyverse)
library(ggplot2)
library(scales)

set.seed(12345)

# Configuration
n_cycles <- 480
n_states <- 19
cycle_length <- 1/12
discount_rate <- 0.035
n_psa <- 1000

# Sample parameters from distributions
sample_parameters <- function() {
  params <- list()
  
  # Progression (Beta)
  params$prog_G2_G3 <- rbeta(1, 60, 940)
  params$prog_G3_G4 <- rbeta(1, 40, 960)
  
  # Exacerbation rates (Gamma)
  params$rate_exac_G2 <- rgamma(1, 25, 166.7)
  params$rate_exac_G3 <- rgamma(1, 36, 120)
  params$rate_exac_G4 <- rgamma(1, 40.3, 73.3)
  
  # CVD rates - no prior CVD (Beta)
  params$rate_MI_G2_noCVD <- rbeta(1, 180, 9820)
  params$rate_MI_G3_noCVD <- rbeta(1, 300, 9700)
  params$rate_MI_G4_noCVD <- rbeta(1, 450, 9550)
  params$rate_stroke_G2_noCVD <- rbeta(1, 120, 9880)
  params$rate_stroke_G3_noCVD <- rbeta(1, 200, 9800)
  params$rate_stroke_G4_noCVD <- rbeta(1, 300, 9700)
  
  # Recurrent CVD (Beta)
  params$rate_MI_postMI <- rbeta(1, 500, 9500)
  params$rate_stroke_postStroke <- rbeta(1, 650, 9350)
  
  # Post-exacerbation risk multipliers (Lognormal)
  params$rr_MI_postExac <- rlnorm(1, log(2.27), 0.096)
  params$rr_stroke_postExac <- rlnorm(1, log(1.83), 0.122)
  params$rr_death_postExac <- rlnorm(1, log(1.50), 0.110)
  
  # Case fatality (Beta)
  params$cfr_MI <- rbeta(1, 850, 7650)
  params$cfr_stroke <- rbeta(1, 1260, 5740)
  params$cfr_exac <- rbeta(1, 250, 4750)
  
  # Mortality (Beta)
  params$mort_COPD_G2 <- rbeta(1, 200, 9800)
  params$mort_COPD_G3 <- rbeta(1, 400, 9600)
  params$mort_COPD_G4 <- rbeta(1, 850, 9150)
  params$mort_background <- rbeta(1, 150, 9850)
  params$mort_excess_postMI <- rbeta(1, 300, 9700)
  params$mort_excess_postStroke <- rbeta(1, 500, 9500)
  
  # Treatment effect (Lognormal)
  params$rr_exac_treatment <- rlnorm(1, log(0.78), 0.062)
  
  # Utilities (Beta)
  params$util_G2_noCVD <- rbeta(1, 156, 44)
  params$util_G3_noCVD <- rbeta(1, 140, 60)
  params$util_G4_noCVD <- rbeta(1, 110, 90)
  params$util_gen_pop <- rbeta(1, 166, 34)
  params$util_postMI_ref <- rbeta(1, 152, 48)
  params$util_postStroke_ref <- rbeta(1, 134, 66)
  params$disutil_exac <- rbeta(1, 26, 174)
  
  # Costs (Gamma)
  params$cost_drug_LAMA <- rgamma(1, 100, 3.5)
  params$cost_drug_LABALAMA <- rgamma(1, 100, 3.08)
  params$cost_maint_G2 <- rgamma(1, 36.8, 1.05)
  params$cost_maint_G3 <- rgamma(1, 36.7, 0.65)
  params$cost_maint_G4 <- rgamma(1, 33.7, 0.32)
  params$cost_postMI_monthly <- rgamma(1, 37.4, 0.53)
  params$cost_postStroke_monthly <- rgamma(1, 23, 0.115)
  params$cost_exac <- rgamma(1, 60.1, 0.028)
  params$cost_MI_acute <- rgamma(1, 51, 0.0066)
  params$cost_stroke_acute <- rgamma(1, 31, 0.0034)
  
  return(params)
}

# Derive combined utilities and cost vectors
calculate_derived_params <- function(params) {
  MI_mult <- params$util_postMI_ref / params$util_gen_pop
  stroke_mult <- params$util_postStroke_ref / params$util_gen_pop
  
  # Combined utilities (multiplicative method)
  params$util_G2_postMI <- params$util_G2_noCVD * MI_mult
  params$util_G2_postStroke <- params$util_G2_noCVD * stroke_mult
  params$util_G3_postMI <- params$util_G3_noCVD * MI_mult
  params$util_G3_postStroke <- params$util_G3_noCVD * stroke_mult
  params$util_G4_postMI <- params$util_G4_noCVD * MI_mult
  params$util_G4_postStroke <- params$util_G4_noCVD * stroke_mult
  
  # Tunnel utilities (additive disutility)
  params$util_G2_noCVD_pE <- max(0, params$util_G2_noCVD - params$disutil_exac)
  params$util_G2_postMI_pE <- max(0, params$util_G2_postMI - params$disutil_exac)
  params$util_G2_postStroke_pE <- max(0, params$util_G2_postStroke - params$disutil_exac)
  params$util_G3_noCVD_pE <- max(0, params$util_G3_noCVD - params$disutil_exac)
  params$util_G3_postMI_pE <- max(0, params$util_G3_postMI - params$disutil_exac)
  params$util_G3_postStroke_pE <- max(0, params$util_G3_postStroke - params$disutil_exac)
  params$util_G4_noCVD_pE <- max(0, params$util_G4_noCVD - params$disutil_exac)
  params$util_G4_postMI_pE <- max(0, params$util_G4_postMI - params$disutil_exac)
  params$util_G4_postStroke_pE <- max(0, params$util_G4_postStroke - params$disutil_exac)
  
  # Cost vectors (18 states + death)
  build_costs <- function(drug) {
    c(params$cost_maint_G2 + drug,
      params$cost_maint_G2 + drug + params$cost_postMI_monthly,
      params$cost_maint_G2 + drug + params$cost_postStroke_monthly,
      params$cost_maint_G2 + drug,
      params$cost_maint_G2 + drug + params$cost_postMI_monthly,
      params$cost_maint_G2 + drug + params$cost_postStroke_monthly,
      params$cost_maint_G3 + drug,
      params$cost_maint_G3 + drug + params$cost_postMI_monthly,
      params$cost_maint_G3 + drug + params$cost_postStroke_monthly,
      params$cost_maint_G3 + drug,
      params$cost_maint_G3 + drug + params$cost_postMI_monthly,
      params$cost_maint_G3 + drug + params$cost_postStroke_monthly,
      params$cost_maint_G4 + drug,
      params$cost_maint_G4 + drug + params$cost_postMI_monthly,
      params$cost_maint_G4 + drug + params$cost_postStroke_monthly,
      params$cost_maint_G4 + drug,
      params$cost_maint_G4 + drug + params$cost_postMI_monthly,
      params$cost_maint_G4 + drug + params$cost_postStroke_monthly,
      0)
  }
  params$state_costs_comp <- build_costs(params$cost_drug_LAMA)
  params$state_costs_intv <- build_costs(params$cost_drug_LABALAMA)
  
  # Utility vector
  params$state_utilities <- c(
    params$util_G2_noCVD, params$util_G2_postMI, params$util_G2_postStroke,
    params$util_G2_noCVD_pE, params$util_G2_postMI_pE, params$util_G2_postStroke_pE,
    params$util_G3_noCVD, params$util_G3_postMI, params$util_G3_postStroke,
    params$util_G3_noCVD_pE, params$util_G3_postMI_pE, params$util_G3_postStroke_pE,
    params$util_G4_noCVD, params$util_G4_postMI, params$util_G4_postStroke,
    params$util_G4_noCVD_pE, params$util_G4_postMI_pE, params$util_G4_postStroke_pE,
    0)
  
  return(params)
}

# Build 19x19 transition matrix
build_transition_matrix <- function(params, arm = "comparator") {
  P <- matrix(0, nrow = n_states, ncol = n_states)
  
  annual_to_monthly <- function(rate) if (rate <= 0) 0 else 1 - exp(-rate / 12)
  exac_mult <- ifelse(arm == "intervention", params$rr_exac_treatment, 1.0)
  
  # State mappings
  tunnel_exit <- c(1,2,3, NA,NA,NA, 7,8,9, NA,NA,NA, 13,14,15, NA,NA,NA, NA)
  stable_to_tunnel <- c(4,5,6, NA,NA,NA, 10,11,12, NA,NA,NA, 16,17,18, NA,NA,NA, NA)
  noCVD_to_postMI <- c(2, NA,NA, 2, NA,NA, 8, NA,NA, 8, NA,NA, 14, NA,NA, 14, NA,NA, NA)
  noCVD_to_postStroke <- c(3, NA,NA, 3, NA,NA, 9, NA,NA, 9, NA,NA, 15, NA,NA, 15, NA,NA, NA)
  prog_target <- c(7,8,9, NA,NA,NA, 13,14,15, NA,NA,NA, NA,NA,NA, NA,NA,NA, NA)
  
  gold_map <- c(rep("II",6), rep("III",6), rep("IV",6), "Death")
  cvd_map <- rep(c("none","postMI","postStroke","none","postMI","postStroke"), 3)
  is_tunnel <- rep(c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE), 3)
  
  for (s in 1:n_states) {
    if (s == 19) { P[19,19] <- 1.0; next }
    
    gold <- gold_map[s]
    cvd <- cvd_map[s]
    tunnel <- is_tunnel[s]
    
    # GOLD-specific rates
    if (gold == "II") {
      rate_exac <- params$rate_exac_G2; rate_MI_base <- params$rate_MI_G2_noCVD
      rate_stroke_base <- params$rate_stroke_G2_noCVD; mort_COPD <- params$mort_COPD_G2
      rate_prog <- params$prog_G2_G3
    } else if (gold == "III") {
      rate_exac <- params$rate_exac_G3; rate_MI_base <- params$rate_MI_G3_noCVD
      rate_stroke_base <- params$rate_stroke_G3_noCVD; mort_COPD <- params$mort_COPD_G3
      rate_prog <- params$prog_G3_G4
    } else {
      rate_exac <- params$rate_exac_G4; rate_MI_base <- params$rate_MI_G4_noCVD
      rate_stroke_base <- params$rate_stroke_G4_noCVD; mort_COPD <- params$mort_COPD_G4
      rate_prog <- 0
    }
    
    # CVD history adjustments
    if (cvd == "postMI") {
      rate_MI <- params$rate_MI_postMI; rate_stroke <- rate_stroke_base * 1.5
      mort_excess <- params$mort_excess_postMI
    } else if (cvd == "postStroke") {
      rate_MI <- rate_MI_base * 1.4; rate_stroke <- params$rate_stroke_postStroke
      mort_excess <- params$mort_excess_postStroke
    } else {
      rate_MI <- rate_MI_base; rate_stroke <- rate_stroke_base; mort_excess <- 0
    }
    
    # Monthly probabilities
    p_MI <- annual_to_monthly(rate_MI)
    p_stroke <- annual_to_monthly(rate_stroke)
    p_exac <- annual_to_monthly(rate_exac) * exac_mult
    p_prog <- annual_to_monthly(rate_prog)
    p_mort <- annual_to_monthly(mort_COPD + params$mort_background + mort_excess)
    
    # Post-exacerbation multipliers
    if (tunnel) {
      p_MI <- p_MI * params$rr_MI_postExac
      p_stroke <- p_stroke * params$rr_stroke_postExac
      p_mort <- p_mort * params$rr_death_postExac
    }
    
    # Case fatality
    p_MI_survive <- p_MI * (1 - params$cfr_MI)
    p_MI_die <- p_MI * params$cfr_MI
    p_stroke_survive <- p_stroke * (1 - params$cfr_stroke)
    p_stroke_die <- p_stroke * params$cfr_stroke
    p_exac_survive <- p_exac * (1 - params$cfr_exac)
    p_exac_die <- p_exac * params$cfr_exac
    
    if (tunnel) {
      P[s, 19] <- min(p_mort + p_MI_die + p_stroke_die, 0.99)
      if (cvd == "none") {
        P[s, noCVD_to_postMI[s]] <- p_MI_survive
        P[s, noCVD_to_postStroke[s]] <- p_stroke_survive
      }
      P[s, tunnel_exit[s]] <- max(1 - sum(P[s,]), 0)
    } else {
      P[s, 19] <- min(p_mort + p_MI_die + p_stroke_die + p_exac_die, 0.99)
      if (cvd == "none" && !is.na(noCVD_to_postMI[s])) {
        P[s, noCVD_to_postMI[s]] <- p_MI_survive
        P[s, noCVD_to_postStroke[s]] <- p_stroke_survive
      }
      if (!is.na(stable_to_tunnel[s])) P[s, stable_to_tunnel[s]] <- p_exac_survive
      if (!is.na(prog_target[s])) P[s, prog_target[s]] <- p_prog
      P[s, s] <- max(1 - sum(P[s,]), 0)
    }
  }
  return(P)
}

# Run Markov trace
run_markov_trace <- function(P, init, n_cycles) {
  trace <- matrix(0, nrow = n_cycles + 1, ncol = n_states)
  trace[1,] <- init
  for (t in 2:(n_cycles + 1)) trace[t,] <- trace[t-1,] %*% P
  return(trace)
}

# Calculate outcomes
calculate_outcomes <- function(trace, params, arm, n_cycles) {
  costs <- if (arm == "comparator") params$state_costs_comp else params$state_costs_intv
  utils <- params$state_utilities
  
  disc <- 1 / (1 + discount_rate) ^ ((0:n_cycles) * cycle_length)
  hc <- c(0.5, rep(1, n_cycles - 1), 0.5)
  
  qaly <- rowSums(trace * utils) * cycle_length
  state_cost <- rowSums(trace * costs)
  
  tunnel_states <- c(4,5,6, 10,11,12, 16,17,18)
  exac_cost <- rowSums(trace[, tunnel_states, drop=FALSE]) * params$cost_exac
  
  cvd_cost <- numeric(n_cycles + 1)
  for (t in 2:(n_cycles + 1)) {
    new_MI <- max(0, sum(trace[t, c(2,8,14)]) - sum(trace[t-1, c(2,8,14)]))
    new_stroke <- max(0, sum(trace[t, c(3,9,15)]) - sum(trace[t-1, c(3,9,15)]))
    cvd_cost[t] <- new_MI * params$cost_MI_acute + new_stroke * params$cost_stroke_acute
  }
  
  total_cost <- state_cost + exac_cost + cvd_cost
  
  return(list(
    qaly_disc = sum(qaly * disc * hc),
    cost_disc = sum(total_cost * disc * hc)
  ))
}

# Initial distribution: 45% GOLD II, 40% GOLD III, 15% GOLD IV; 80% no CVD
init <- rep(0, n_states)
init[1] <- 0.45 * 0.80; init[2] <- 0.45 * 0.12; init[3] <- 0.45 * 0.08
init[7] <- 0.40 * 0.80; init[8] <- 0.40 * 0.12; init[9] <- 0.40 * 0.08
init[13] <- 0.15 * 0.80; init[14] <- 0.15 * 0.12; init[15] <- 0.15 * 0.08

# Run PSA
cat("Running PSA (n =", n_psa, ")...\n")
psa_results <- data.frame(
  iteration = 1:n_psa, cost_comp = NA, cost_intv = NA,
  qaly_comp = NA, qaly_intv = NA, incr_cost = NA, incr_qaly = NA
)

pb <- txtProgressBar(min = 0, max = n_psa, style = 3)
for (i in 1:n_psa) {
  params <- calculate_derived_params(sample_parameters())
  
  P_comp <- build_transition_matrix(params, "comparator")
  P_intv <- build_transition_matrix(params, "intervention")
  
  trace_comp <- run_markov_trace(P_comp, init, n_cycles)
  trace_intv <- run_markov_trace(P_intv, init, n_cycles)
  
  out_comp <- calculate_outcomes(trace_comp, params, "comparator", n_cycles)
  out_intv <- calculate_outcomes(trace_intv, params, "intervention", n_cycles)
  
  psa_results$cost_comp[i] <- out_comp$cost_disc
  psa_results$cost_intv[i] <- out_intv$cost_disc
  psa_results$qaly_comp[i] <- out_comp$qaly_disc
  psa_results$qaly_intv[i] <- out_intv$qaly_disc
  psa_results$incr_cost[i] <- out_intv$cost_disc - out_comp$cost_disc
  psa_results$incr_qaly[i] <- out_intv$qaly_disc - out_comp$qaly_disc
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Summary
cat("\n\nPSA RESULTS\n")
cat(sprintf("Incr QALYs: %.4f (95%% CI: %.4f to %.4f)\n",
            mean(psa_results$incr_qaly),
            quantile(psa_results$incr_qaly, 0.025),
            quantile(psa_results$incr_qaly, 0.975)))
cat(sprintf("Incr Costs: £%.0f (95%% CI: £%.0f to £%.0f)\n",
            mean(psa_results$incr_cost),
            quantile(psa_results$incr_cost, 0.025),
            quantile(psa_results$incr_cost, 0.975)))

prob_dom <- mean(psa_results$incr_qaly > 0 & psa_results$incr_cost < 0)
prob_ce_20k <- mean(psa_results$incr_qaly > 0 & 
                      (psa_results$incr_cost < 0 | psa_results$incr_cost/psa_results$incr_qaly < 20000))
prob_ce_30k <- mean(psa_results$incr_qaly > 0 & 
                      (psa_results$incr_cost < 0 | psa_results$incr_cost/psa_results$incr_qaly < 30000))

cat(sprintf("P(dominant): %.1f%%\n", prob_dom * 100))
cat(sprintf("P(CE at £20k): %.1f%%\n", prob_ce_20k * 100))
cat(sprintf("P(CE at £30k): %.1f%%\n", prob_ce_30k * 100))

# CE Plane
ce_plane <- ggplot(psa_results, aes(x = incr_qaly, y = incr_cost)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 20000, linetype = "dashed", color = "darkgreen") +
  geom_abline(intercept = 0, slope = 30000, linetype = "dashed", color = "orange") +
  labs(title = "Cost-Effectiveness Plane", x = "Incremental QALYs", y = "Incremental Costs (£)") +
  scale_y_continuous(labels = comma) +
  theme_minimal()

ggsave("/mnt/user-data/outputs/CE_plane.png", ce_plane, width = 8, height = 6, dpi = 150)

# CEAC
wtp_range <- seq(0, 50000, 1000)
prob_ce <- sapply(wtp_range, function(wtp) {
  mean(psa_results$incr_qaly > 0 & 
         (psa_results$incr_cost < 0 | psa_results$incr_cost/psa_results$incr_qaly < wtp))
})

ceac <- ggplot(data.frame(wtp = wtp_range, prob = prob_ce), aes(x = wtp, y = prob)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_vline(xintercept = 20000, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = 30000, linetype = "dashed", color = "orange") +
  labs(title = "Cost-Effectiveness Acceptability Curve", 
       x = "WTP Threshold (£/QALY)", y = "Probability Cost-Effective") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(limits = c(0, 1), labels = percent) +
  theme_minimal()

ggsave("/mnt/user-data/outputs/CEAC.png", ceac, width = 8, height = 6, dpi = 150)

# Save results
write.csv(psa_results, "/mnt/user-data/outputs/PSA_results.csv", row.names = FALSE)
cat("\nOutputs saved to /mnt/user-data/outputs/\n")