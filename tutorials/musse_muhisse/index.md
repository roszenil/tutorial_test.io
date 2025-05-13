---
title: Hypothesis testing with the Multiple State-dependent Speciation and Extinction (MuSSE) and Multiple Hidden State-dependent Speciation and Extinction (MuHiSSE) models in RevBayes
subtitle: 
authors:  Jenna McCullough
level: 0
order: 2
index: true
redirect: false
include_files:
    - musse_muhisse/data/combined_traitdatabase_0123.tsv
    - musse_muhisse/data/MCC_AFRICANRESTIOS.tre
    - tutorial_structure/scripts/test.Rev
    - tutorial_structure/scripts/test.Rmd



---

{% section Overview %}
The purpose of this tutorial is to give a detailed description of how to use Multiple State-dependent Speciation and Extinction (MuSSE) and Multiple Hidden State-dependent Speciation and Extinction (MuHiSSE) models to test specific hypotheses related to state-dependent diversification in RevBayes. Because the goal is to describe the critical thinking steps that goes into hypothesis testing with Bayesian SSEs, this tutorial also describes a case study for Restios (described in Zenil-Ferguson et al. XXXX) and the detailed interpretation of the output files. 

{% subsection MuSSE and MuHiSSE %}

MuSSE can accomodate three or more states in which to test hypotheses about state-dependent diversification. MuHiSSE is a model extension of the MuSSE model by its addition of a hidden state. These models incorporate more than two states; we show four states below. Each state has its own speciation (λ, lambda) and extinction (μ, mu) parameters as well as anagenetic transitions between states (q parameters). MuHiSSE additionally has parameters that govern the transitions between the A and B states (alpha and beta). It is important to mention here that, similar to Binary Speciation and Extinction (BiSSE) and  Hidden State Speciation and Extinction (HiSSE), these models only allow for a trait change to occur along a branch (anagenetic trait evolution, denoted by the q transition parameters) rather than at a speciation event (cladogenetic trait evolution). If this is something you would like to implement in your own work, see the Cladogenetic Speciation and Exinction (ClaSSE, Goldberg and Igic 2012).

{% figure example %}
<img src="figures/02_MuSSE-MuHiSSE.pdf" width="600">
{% figcaption %}
Differences in the number of states and parameters between MuSSE and MuHiSSE models. Each state has a speciation (λ, lambda), extinction (μ, mu), anagenetic transitions between states (q) and hidden states (alpha and beta) rates. 
{% endfigcaption %}
{% endfigure %}

In RevBayes, it is simple to control the transition matrix if your own empirical system required it. For example, if it was biologically unrealistic for an anagenetic transition to occur from state 0 to state 2, it would be possible to set the transition matrix to only allow transitions into state two from either state 1 or 3. This is discussed more below.

In this tutorial, we will not be going through the step by step guide for MuSSE because taking the time to perform their more complex, hidden state extensions in a Bayesian framework inherently includes their more simpler versions. In other words, we can interpret results from only the A states (0A, 1A, 2A, and 3A) from the MuHiSSE if we wanted to know the results of a non-hidden state model. But by running the more complex model in RevBayes, we have the information for both and the confidence that our results will be less likely to be hampered by false positive inferences. 

The rev code between MuSSE and MuHiSSE is not significanty different but a MuSSE-specific script `musse.rev` is be available in at the top of this webpage. 

{% subsection Tensorphylo %}

Because Bayesian SSE models are computationally intensive, especially for large trees, TensorPhylo is a program that can improve computational time. Because it might be impractical for some users to implement TensorPhylo, we have versions of this tutorial with and without Tensorphylo.

NOTE: should we just include the bit about Tensorphylo at the bottom of the tutorial instead? 

{% subsection Testing Hypotheses and Overall Workflow %}

SSE models are powerful methods to explore the relationship of traits to the tempo and mode of evolution. To use them most effectively, we need to have specific hypotheses that can only come from understanding the biology of your study system. 

The data we will be analyzing for this tutorial is the same as discussed in Zenil-Ferguson et al. (XXXXXXX): a reanalyzed dataset of restios (Polaes:Restionaceae from Bouchenak-Khelladi and Linder 2017), which are a family of bamboo-like flowering plants mostly found in the Cape Floristic Region of South Africa. Bouchenak-Khelladi and Linder (2017) sought to understand how two binary habtiat traits, precipitation (dry or wetlands) and elevation (costal or mountainous), contribute to faster rates of diversification. They also combined these two traits into four states for a MuSSE model (coastal dryland, coastal wetlands, mountainous drylands, and mountainous wetlands), which is what we will be running in this tutorial. 

They had two hypotheses about what contributed to higher diversification:

1. The older mountainous wetland habitats are more fragmented and have promoted lineage diversification because of higher rates of turnover.

2. The recent coastal dryland habitats have had much more connectivity and geographic area in which to diversify.

The authors fitted BiSSE, HiSSE, ClaSSE, and MuSSE in the R packages hisse and diversitree, both likelihood methods, and determined the best model with the Akaike Information Criterion (AIC; See the likelihood workflow below). When considering elevation and precipitation independently (BiSSE, HiSSE, and ClaSSE), Bouchenak-Khelladi and Linder (2017) found no relationship with net diversification for either trait. When they combined the traits in a 4-state MuSSE model, the best-fit model was the one that had the same speciation but different extinction rates. Overall, Bouchenak-Khelladi and Linder (2017) concluded that no habitat category had higher net diversification rates, thus not supporting their original hypotheses.

{% figure example %}
<img src="figures/likelihood-bayesian-comparison.pdf" width="600">
{% figcaption %}
This is a simplified workflow for implementing Likelihood vs. Bayesian SSE model 
{% endfigcaption %}
{% endfigure %}

NOTE TO ROSANA ET AL.: Is this figure worth including here? I thought it would be good to again reiterate how the workflow of Bayesian SSEs is different than Likelihood. I'm not going to edit the Bayesian workflow side to match the number of states for the MuHiSSE model yet in case people do not think it is necessary. 

Instead of running a MuHiSSE model in which we set expectation for rates (i.e., the difference in the models A–D in the workflow figure above), we will analyze one model, a MuHiSSE allowing with free rates. We will then in

{% subsection Data preparation %}

Begin by creating a project _directory_ (also referred to as a folder) on your computer, titled `musse` with two subdirectories inside, one called _data_ and the other called _scripts_. See the box called `Data files and scripts` in the upper left-hand corner of this webpage and download the files. Put `combined_traitdatabase_0123.tsv` and `MCC_AFRICANRESTIOS.tre` into your `data` directory and `test.Rev` and `test.Rmd` into your `scripts` directory.

{% subsection Setup %}
This tutorial uses RevBayes v1.2.5. Earlier versions may not run properly. 

The two necessary data files that MuHiSSE 



To begin, load in the sequences using the ‘readDiscreteCharacterData()‘
function.
```
data_cox2 = readDiscreteCharacterData("data/primates_and_galeopterus_cox2.nex")
data_cytb = readDiscreteCharacterData("data/primates_and_galeopterus_cytb.nex")
```

