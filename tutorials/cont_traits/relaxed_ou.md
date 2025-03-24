---
title: Relaxed Ornstein-Uhlenbeck Models
subtitle: Estimating lineage-specific optima under Ornstein-Uhlenbeck evolution
authors: Michael R. May and Sebastian Höhna
level: 6
order: 1.8
prerequisites:
- intro
- intro/revgadgets
- mcmc
- cont_traits/cont_trait_intro
- cont_traits/simple_ou
index: true
redirect: false
include_all: false
include_files:
- data/primates_tree.nex
- data/primates_cont_traits.nex
- scripts/mcmc_relaxed_OU.Rev
- scripts/plot_relaxed_OU.R
---

{% section Estimating Branch-Specific Evolutionary Optima %}

This tutorial demonstrates how to specify an Ornstein-Uhlenbeck model where the optimal phenotype is allowed to vary over branches of a time-calibrated phylogeny {% cite Uyeda2014 %}. We provide the probabilistic graphical model representation of each component for this tutorial. After specifying the model, you will estimate the parameters of branch-specific Ornstein-Uhlenbeck evolution using reversible-jump Markov chain Monte Carlo (rjMCMC).

{% include_relative modules/relaxed_OU_parameter_estimation.md %}

{% include_relative modules/relaxed_OU_exercise_1.md %}
