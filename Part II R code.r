# Key Multiplicity Issues in Clinical Trials (Part II)
# Online training course by Alex Dmitrienko
# Downloaded from https://github.com/medianasoft/MultCompSAS

library(Mediana)
library(gsDesign)

##############################################################

# Module C

# Computation of adjusted p-values for the Hochberg-based 
# parallel gatekeeping procedure

##############################################################

# Example 2: Type 2 diabetes mellitus trial

rawp=c(0.0082,0.0174,0.0202) 
families=families(family1=c(1, 2),
                  family2=c(3))
component.procedures=
              families(family1 ="HochbergAdj",
                       family2="HochbergAdj")
gamma=families(family1=0.7,
                 family2=1)
adjp=AdjustPvalues(rawp,
          proc="ParallelGatekeepingAdj",
          par=parameters(family=families,
                         proc=component.procedures,
                         gamma=gamma))
round(adjp, 4)

##############################################################

# Module D

# Computation of adjusted p-values for the Hochberg-based 
# general gatekeeping procedure

##############################################################

# Example 3: Schizophrenia trial

rawp=c(0.0101,0.0233,0.0022,0.0167) 
families=families(family1=c(1, 2),
                  family2=c(3, 4))
component.procedures=
              families(family1 ="HochbergAdj",
                       family2="HochbergAdj")
gamma=families(family1=0.7,
                 family2=1)
adjp=AdjustPvalues(rawp,
          proc="MultipleSequenceGatekeepingAdj",
          par=parameters(family=families,
                         proc=component.procedures,
                         gamma=gamma))
round(adjp, 4)

##############################################################

# Module E

# Computation of critical values in a group-sequential design

##############################################################

# Example 4: Prostate cancer trial

gsDesign(k=2, timing=c(0.5, 1), test.type=1, 
       sfu="OF", alpha=0.025, beta=0.1)

gsDesign(k=2, timing=c(0.5, 1), test.type=1, 
       sfu="OF", alpha=0.025 * 0.8, beta=0.1)

gsDesign(k=2, timing=c(0.5, 1), test.type=1, 
       sfu="OF", alpha=0.025 * 0.8, beta=0.1)

##############################################################

# Module G

# Power calculations

##############################################################

# Example 3: Schizophrenia trial

# Outcome distribution parameters
placebo.par = parameters(
              parameters(mean = -12, sd = 20),
              parameters(mean = -0.8, sd = 1))

dosel.par = parameters(
              parameters(mean = -18, sd = 20),
              parameters(mean = -1.2, sd = 1))

doseh.par = parameters(
              parameters(mean = -20, sd = 20),
              parameters(mean = -1.3, sd = 1))

# Correlation between two endpoints
corr.matrix = matrix(c(1.0, 0.3,
                       0.3, 1.0), 2, 2)

# Variable types
var.type = list("NormalDist", "NormalDist")

# Outcome parameters
outcome.placebo = parameters(type = var.type,
                              par = placebo.par,
                              corr = corr.matrix)
outcome.dosel = parameters(type = var.type,
                            par = dosel.par,
                            corr = corr.matrix)
outcome.doseh = parameters(type = var.type,
                            par = doseh.par,
                            corr = corr.matrix)

# Data model
ex3.data.model = DataModel() +
  OutcomeDist(outcome.dist = "MVNormalDist") +
  SampleSize(c(110, 140)) +
  Sample(id = list("Placebo (P)", "Placebo (S)"),
         outcome.par = parameters(outcome.placebo)) +
  Sample(id = list("Dose L (P)", "Dose L (S)"),
         outcome.par = parameters(outcome.dosel)) +
  Sample(id = list("Dose H (P)", "Dose H (S)"),
         outcome.par = parameters(outcome.doseh))

# Parameters of the gatekeeping procedure procedure (multiple-sequence gatekeeping procedure)

# Families of hypotheses
families = families(family1 = c(1, 2),
                  family2 = c(3, 4))

# Component procedures for each family
component.procedures = 
        families(family1 ="HochbergAdj",
                 family2 = "HochbergAdj")

# Truncation parameter for each family
gamma = families(family1 = 0.7,
                 family2 = 1)

# Analysis model
ex3.analysis.model = AnalysisModel() +
  MultAdjProc(proc = "MultipleSequenceGatekeepingAdj",
              par = parameters(family = families,
                               proc = component.procedures,
                               gamma = gamma),
              tests = tests("Placebo vs Dose L (P)",
                  "Placebo vs Dose H (P)",
                  "Placebo vs Dose L (S)",
                  "Placebo vs Dose H (S)")) +
  Test(id = "Placebo vs Dose L (P)",
       method = "TTest",
       samples = samples("Dose L (P)", "Placebo (P)")) +
  Test(id = "Placebo vs Dose H (P)",
       method = "TTest",
       samples = samples("Dose H (P)", "Placebo (P)")) +
  Test(id = "Placebo vs Dose L (S)",
       method = "TTest",
       samples = samples("Dose L (S)", "Placebo (S)")) +
  Test(id = "Placebo vs Dose H (S)",
       method = "TTest",
       samples = samples("Dose L (S)", "Placebo (S)"))

# Evaluation model
ex3.evaluation.model = EvaluationModel() +
  Criterion(id = "Disjunctive power (P)",
            method = "DisjunctivePower",
            tests = tests("Placebo vs Dose L (P)",
                          "Placebo vs Dose H (P)"),
            labels = tests("Disjunctive power"),
            par = parameters(alpha = 0.025)) +
  Criterion(id = "Disjunctive power (S)",
            method = "DisjunctivePower",
            tests = tests("Placebo vs Dose L (S)",
                          "Placebo vs Dose H (S)"),
            labels = tests("Disjunctive power"),
            par = parameters(alpha = 0.025))

# Simulation Parameters
ex3.sim.parameters =  
         SimParameters(n.sims = 10000,
                       proc.load = "full",
                       seed = 42938001)

# Perform clinical scenario evaluation
ex3.results = CSE(ex3.data.model,
                          ex3.analysis.model,
                          ex3.evaluation.model,
                          ex3.sim.parameters)

summary(ex3.results)
