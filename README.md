# Thesis Course 2025/2026

## A Unified Statistical Framework for Efficacy-Integrated Dose Optimisation Designs in Early-Phase Oncology Trials

This repository contains scripts for generating simulated trial data, along with the corresponding datasets and analysis results, to compare efficacy-integrated dose-finding designs in early-phase oncology trials.

## Project overview

The project develops and evaluates a unified simulation framework for comparing model-assisted dose-finding methods that jointly consider toxicity and efficacy. The focus is on how different designs allocate patients during the trial and how accurately they identify the optimal biological dose (OBD) at the end of the trial.

## Designs considered

The following designs are included in this repository:

- **TEPI** - Toxicity and Efficacy Probability Interval
- **STEIN** - Simple Toxicity and Efficacy Interval
- **BOIN12** - Bayesian Optimal Interval Phase I/II
- **RWR** - Random Walk Rule comparator

## Outcomes

Patient outcomes are represented using a three-outcome framework:

- **Neutral**: no toxicity and no efficacy
- **Success**: efficacy without toxicity
- **Toxicity**: toxicity, regardless of efficacy

## Data-generating framework

The simulation study is based on benchmark scenarios motivated by the **LORDs** framework for trinomial dose-finding outcomes. The following data-generating settings are considered:

- **Base case**: logit/logit
- **Misspecification 1**: probit/probit
- **Misspecification 2**: quadratic efficacy with logistic toxicity

## Repository structure

```text

├── SRC/
│   ├── simulation-tepi_stein_rwr.R
│   ├── simulation-boin12.R
│   ├── plots_thesis.R
├── Results/
│   ├── tables/
│   └── figures/
|   └── Simulated Data/
└── README.Rmd

```
