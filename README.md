# p53-Mdm2 ODE Model for Cisplatin Sensitivity Prediction in NSCLC

## Overview

This project investigates whether a mechanistic ODE model of p53-Mdm2 
dynamics can predict cisplatin sensitivity in NSCLC cell lines, and how 
it compares against ML approaches trained on the same data.

## Biological Background

p53 is a tumour suppressor protein regulated by Mdm2 through a negative 
feedback loop. During DNA damage, phosphorylation of both proteins 
disrupts Mdm2-mediated degradation, causing p53 accumulation and 
triggering apoptosis. TP53 is commonly mutated in NSCLC, and drugs like 
cisplatin work by inducing p53-dependent apoptosis — making this pathway 
a key factor in treatment response.

## Research Questions

1. Can an ODE model of p53-Mdm2 dynamics predict cisplatin sensitivity 
in NSCLC cell line data?
2. How do ODE model predictions compare against ML approaches trained 
on the same genomic features?

## Data

- **GDSC2** — cisplatin IC50 values for NSCLC cell lines
- **DepMap** — gene expression and somatic mutation data 
  (OmicsExpressionProteinCodingGenesTPMLogp1.csv, 
  OmicsSomaticMutationsMatrixDamaging.csv)
- 57 LUAD/LUSC cell lines retained after filtering and matching

## Methods

### ODE Model
- p53-Mdm2 interaction modelled as a system of ODEs using the deSolve 
  package in R
- Fixed DNA damage signal DS = 0.8
- Each cell line simulated for 72 hours
- Four features extracted: peak p53, mean p53, apoptosis score, 
  late phase p53
- Ridge regression used to predict IC50 from ODE features

### ML Models
Three models trained on genomic features (TP53 mutation, Mdm2 
expression, ATM expression and mutation status):
- Elastic Net
- Random Forest
- Neural Network

## Results

| Model | Pearson r |
|-------|-----------|
| ODE + Ridge | 0.309 |
| Elastic Net | 0.321 |
| Neural Network | 0.560 |
| Random Forest | 0.669 |

Random Forest outperformed all models. The ODE model performed 
comparably to Elastic Net, suggesting the mechanistic model captures 
similar linear information from the four genomic features despite using 
no direct genomic input — only simulated p53 pathway dynamics.

The ODE model correctly separates TP53 wild-type and mutant cell lines, 
confirming that p53 pathway dynamics encode biologically meaningful 
signal for cisplatin response.

## Key Finding

A mechanistic model built from first-principles biology achieves 
comparable predictive performance to a linear ML model, while offering 
interpretability that data-driven approaches cannot provide. The two 
approaches are complementary rather than competitive.

## Tech Stack

R, deSolve, glmnet, randomForest, nnet, ggplot2

## Reference

Assignment II — MEIN40330 AI for Personalised Medicine, University 
College Dublin
