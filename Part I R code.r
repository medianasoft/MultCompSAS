# Key Multiplicity Issues in Clinical Trials (Part I)
# Online training course by Alex Dmitrienko
# Downloaded from https://github.com/medianasoft/MultCompSAS

library(Mediana)

##############################################################

# Module C

# Computation of adjusted p-values for stepwise procedures 
# with a data-driven testing sequence

##############################################################

# Example 4: Type 2 diabetes trial

rawp = c(0.0111,0.0065,0.0293)
AdjustPvalues(rawp, proc = "HolmAdj")

# Example 1: Prostate cancer trial

rawp = c(0.0102,0.0181)
weight = c(0.8,0.2)
AdjustPvalues(rawp, proc = "HolmAdj",
              par = parameters(weight = weight))

##############################################################

# Module D

# Computation of adjusted p-values for stepwise procedures 
# with a pre-specified testing sequence

##############################################################

# Example 4: Type 2 diabetes trial

# Fallback procedure 

rawp=c(0.0061,0.0233,0.0098)
weight=c(1/3,1/3,1/3)
transition = matrix(c(0, 1, 0,
                      0, 0, 1,
                      0, 0, 0), 3, 3, byrow = TRUE)

AdjustPvalues(rawp, proc = "ChainAdj",
  par = parameters(weight = weight,
  transition = transition))

# Chain procedure 

rawp=c(0.0061,0.0233,0.0098)
weight=c(1/3,1/3,1/3)
transition = matrix(c(0, 0.5, 0.5,
                      0, 0, 1,
                      0, 0, 0), 3, 3, byrow = TRUE)

AdjustPvalues(rawp, proc = "ChainAdj",
  par = parameters(weight = weight,
  transition = transition))

##############################################################

# Module E

# Computation of adjusted p-values for parametric procedures

##############################################################

# Example 4: Type 2 diabetes trial

stat=c(2.64, 1.93, 2.31)
n=90
p=1-pt(stat, df=2*(n-1))
adjp=AdjustPvalues(p, proc="DunnettAdj",
  par = parameters(n=n))
round(adjp, 4) 

##############################################################

# Module F

# Computation of lower limits of one-sided simultaneous confidence 
# intervals for nonparametric and parametric procedures

##############################################################

# Example 4: Type 2 diabetes trial

est=c(0.63, 0.46, 0.55)

ci=AdjustCIs(est, proc="HolmAdj",
  par = parameters(sd=rep(1.6, 3), n=90, covprob=0.975))
round(ci, 3) 

##############################################################

# Module G

# Power calculations in clinical trials with multiple objectives

##############################################################

# Example 4: Type 2 diabetes trial

# Analytical approach to sample size and power calculations

library(mvtnorm)

alpha = 0.025
beta = 0.2
z_alpha = qnorm(1 - alpha)
z_beta = qnorm(1 - beta)
theta = 0.3

# Standard sample size formula
n_standard = ceiling(2 * (z_alpha + z_beta)^2 / theta^2)
n_standard

# Modified sample size formula
z_alpha = qnorm(1 - alpha / 3)
n_modified = ceiling(2 * (z_alpha + z_beta)^2 / theta^2)
n_modified

# Analytical power calculations without adjustment

# Modified sample size formula
n = n_modified

z_alpha = qnorm(1 - alpha)

# Marginal power of each test
pnorm(sqrt(n / 2) * theta - z_alpha) 

# Disjunctive power
corr = matrix(c(1, 0.5, 0.5,
                0.5, 1, 0.5,
                0.5, 0.5, 1), 3, 3)
1 - pmvnorm(lower = -Inf, upper = rep(z_alpha, 3), mean = rep(sqrt(n / 2) * theta, 3), corr = corr)

# Analytical power calculations with a naive Bonferroni adjustment

z_alpha = qnorm(1 - alpha / 3)

# Marginal power of each test
pnorm(sqrt(n / 2) * theta - z_alpha) 

# Disjunctive power
corr = matrix(c(1, 0.5, 0.5,
                0.5, 1, 0.5,
                0.5, 0.5, 1), 3, 3)
1 - pmvnorm(lower = -Inf, upper = rep(z_alpha, 3), mean = rep(sqrt(n / 2) * theta, 3), corr = corr)

# Example 4: Type 2 diabetes trial

# Mediana code

# Data model
outcome.placebo = parameters(mean = 0, sd = 1)
outcome.dose1 = parameters(mean = 0.3, sd = 1)
outcome.dose2 = parameters(mean = 0.3, sd = 1)
outcome.dose3 = parameters(mean = 0.3, sd = 1)

ex4.data.model = DataModel() +
  OutcomeDist(outcome.dist = "NormalDist") +
  SampleSize(seq(130, 150, 5)) +
  Sample(id = "Placebo",
    outcome.par = parameters(outcome.placebo)) +
  Sample(id = "Dose 1",
    outcome.par = parameters(outcome.dose1)) +
  Sample(id = "Dose 2",
    outcome.par = parameters(outcome.dose2)) +
  Sample(id = "Dose 3",
    outcome.par = parameters(outcome.dose3))

# Analysis model
ex4.analysis.model = AnalysisModel() +
  MultAdjProc(proc = "BonferroniAdj") +
  Test(id = "Placebo vs Dose 1",
    samples = samples("Placebo", "Dose 1"),
    method = "TTest") +
  Test(id = "Placebo vs Dose 2",
    samples = samples ("Placebo", "Dose 2"),
    method = "TTest") +
  Test(id = "Placebo vs Dose 3",
    samples = samples("Placebo", "Dose 3"),
    method = "TTest")

# Evaluation model
ex4.evaluation.model = EvaluationModel() +
  Criterion(id = "Marginal power",
    method = "MarginalPower",
    tests = tests("Placebo vs Dose 1",
                  "Placebo vs Dose 2",
                  "Placebo vs Dose 3"),
    labels = c("Placebo vs Dose 1",
               "Placebo vs Dose 2",
               "Placebo vs Dose 3"),
    par = parameters(alpha = 0.025)) +
  Criterion(id = "Disjunctive power",
    method = "DisjunctivePower",
    tests = tests("Placebo vs Dose 1",
                  "Placebo vs Dose 2",
                  "Placebo vs Dose 3"),
    labels = "Disjunctive power",
    par = parameters(alpha = 0.025))  

# Simulation parameters
ex4.sim.parameters = 
    SimParameters(n.sims = 10000,
                  proc.load = "full",
                  seed = 42938001)

# Perform clinical scenario evaluation
ex4.results = CSE(ex4.data.model,
                  ex4.analysis.model,
                  ex4.evaluation.model,
                  ex4.sim.parameters)

# Simple summary of the simulation results
summary(ex4.results)


# Example 5: Non-small-cell lung cancer trial

# Mediana code

# Data model
outcome.pn = parameters(rate = log(2) / 11)
outcome.pp = parameters(rate = log(2) / 11)
outcome.tn = parameters(rate = log(2) / 12.5)
outcome.tp = parameters(rate = log(2) / 15)

sample.size.total = c(530)
sample.size.pn = as.list(0.2 * sample.size.total)
sample.size.pp = as.list(0.3 * sample.size.total)
sample.size.tn = as.list(0.2 * sample.size.total)
sample.size.tp = as.list(0.3 * sample.size.total)

ex5.data.model = DataModel() +
  OutcomeDist(outcome.dist = "ExpoDist") +
  Sample(id = "Placebo (marker-negative)",
         sample.size = sample.size.pn,
         outcome.par = parameters(outcome.pn)) +
  Sample(id = "Placebo (marker-positive)",
         sample.size = sample.size.pp,
         outcome.par = parameters(outcome.pp)) +
  Sample(id = "Treatment (marker-negative)",
         sample.size = sample.size.tn,
         outcome.par = parameters(outcome.tn)) +
  Sample(id = "Treatment (marker-positive)",
         sample.size = sample.size.tp,
         outcome.par = parameters(outcome.tp))

# Analysis model
ex5.analysis.model = AnalysisModel() +
  MultAdjProc(proc = "HochbergAdj") +
  Test(id = "Overall population test",
       samples = samples(c("Placebo (marker-negative)", 
                           "Placebo (marker-positive)"),
                         c("Treatment (marker-negative)", 
                           "Treatment (marker-positive)")),
       method = "LogrankTest") +
  Test(id = "Marker-positive subpopulation test",
       samples = samples("Placebo (marker-positive)", 
                         "Treatment (marker-positive)"),
       method = "LogrankTest")

# Evaluation model
ex5.evaluation.model = EvaluationModel() +
  Criterion(id = "Marginal power",
            method = "MarginalPower",
            tests = tests("Overall population test",
                          "Marker-positive subpopulation test"),
            labels = c("Overall population test",
                       "Marker-positive subpopulation test"),
            par = parameters(alpha = 0.025)) +
  Criterion(id = "Disjunctive power",
            method = "DisjunctivePower",
            tests = tests("Overall population test",
                          "Marker-positive subpopulation test"),
            labels = "Disjunctive power",
            par = parameters(alpha = 0.025)) 

# Simulation parameters
ex5.sim.parameters = 
    SimParameters(n.sims = 10000,
                  proc.load = "full",
                  seed = 42938001)

# Perform clinical scenario evaluation
ex5.results = CSE(ex5.data.model,
                  ex5.analysis.model,
                  ex5.evaluation.model,
                  ex5.sim.parameters)

# Simple summary of the simulation results
summary(ex5.results)
