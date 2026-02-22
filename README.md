# COPD_economic_model
When health technology assessment bodies such as NICE evaluate new treatments for chronic obstructive pulmonary disease (COPD), they typically commission economic models that simulate patient trajectories through defined health states. These models track disease progression from mild to moderate to severe disease, capture exacerbation events, account for mortality, and ultimately calculate costs and quality-adjusted life years (QALYs) for the intervention versus standard care. The standard architecture for such models includes COPD severity defined by GOLD stage, exacerbation events modelled as transient episodes with associated costs and utility decrements, background mortality adjusted for disease severity, and treatment effects derived from randomised controlled trials.
This framework has guided resource allocation decisions for COPD therapies for over two decades. However, a growing body of epidemiological evidence suggests that these models may systematically underestimate the value of therapies that reduce exacerbations by failing to capture an important downstream consequence: the transient elevation in cardiovascular risk that follows acute respiratory events.
This analysis presents a 19-state Markov model designed to capture what standard approaches miss. By incorporating the dynamic interaction between COPD exacerbations and cardiovascular events, the model quantifies the additional value of exacerbation-reducing therapies that current health technology assessment submissions may not always recognise.
This model captures this interaction through:
- **19 health states** defined by GOLD stage (II-IV), CVD history (none/post-MI/post-stroke), 
  and post-exacerbation status (stable/tunnel)
- **Tunnel states** representing the 30-day post-exacerbation period with elevated 
  cardiovascular transition probabilities
- **Indirect treatment benefit**: exacerbation reduction → less tunnel time → fewer CVD events

## Key Findings

| Outcome | Value |
|---------|-------|
| Incremental QALYs | +0.114 |
| Incremental costs | -£174 |
| ICER | Dominant |
| Probability cost-effective (£20k threshold) | 99.8% |
| CVD contribution to QALY gain | 38% |
| CVD contribution to cost savings | 81% |

## Repository Contents

- `COPD_CVD_Model_19State_Full.xlsx` - Excel workbook with transition matrices and Markov trace
- `COPD_CVD_Parameters.xlsx` - Parameter values with sources and distribution specifications
- `OWSA_simplified.py` - One-way sensitivity analysis with tornado diagram
- `COPD_CVD_PSA_Model.R` - Probabilistic sensitivity analysis (1,000 iterations)
- `/outputs` - OWSA tornado diagram, CE plane, CEAC, results tables

## Model Specifications

| Specification | Value |
|---------------|-------|
| Perspective | NHS and Personal Social Services |
| Time horizon | Lifetime (40 years) |
| Cycle length | 1 month |
| Discount rate | 3.5% (costs and outcomes) |
| Half-cycle correction | Applied |

## Topics/Tags for GitHub
health-economics, cost-effectiveness, markov-model, copd, cardiovascular, 
health-technology-assessment, qaly, icer, sensitivity-analysis, nice
```
health-economics, cost-effectiveness, markov-model, copd, cardiovascular, 
health-technology-assessment, qaly, icer, sensitivity-analysis, nice
