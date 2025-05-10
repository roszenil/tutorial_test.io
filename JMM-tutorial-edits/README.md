# Bayesian Hypothesis Testing with SSE models in RevBayes 
## Tutorial overview
The purpose of this tutorial is to describe Bayesian hypothesis testing with SSE models, specifically with Hidden State-dependent Speciation and Extinction (HiSSE) and Hidden Multiple State-dependent Speciation and Extinction (MuHiSSE) models in RevBayes. 

Both HiSSE and MuHiSSE estimate differences in state-dependent diversification but differ in the number of allowed states: HiSSE requires binary states (0 and 1) and MuHiSSE allows for more than two states. HiSSE and MuHiSSE are the hidden state extensions, respectively, of the Binary State-dependent Speciation and Extinction (BiSSE) and Multiple State-dependent Speciation and Extinction (MuSSE) models. It is important to mention here that all these models only allow for a trait change to occur along a branch (anagenetic trait evolution, denoted by the q transition parameters) rather than at a speciation event (cladogenetic trait evolution). If this is something you would like to implement in your own work, see the Cladogenetic State-dependent Speciation and Exinction (ClaSSE, Goldberg and Igic 2012). 

<img width="975" alt="Screen Shot 2025-05-07 at 2 32 37 PM" src="https://github.com/user-attachments/assets/eb26d248-be98-42b3-9746-948415ea3385" />
HiSSE is a model extension of the BiSSE model by its addition of a hidden state. Both have two states, 0 and 1, and HiSSE has the additional layer of complexity by having A and B versions of the two states. Each state has its own speciation (λ, lambda) and extinction (μ, mu) parameters as well as transitions between states (q parameters). HiSSE additionally has parameters that govern the transitions between the A and B states (alpha and beta). 
<img width="1184" alt="Screen Shot 2025-05-07 at 2 58 43 PM" src="https://github.com/user-attachments/assets/256a7b8d-43be-4a4d-b547-4698d0c696a5" />


MuHiSSE is a model extension of the MuSSE model by its addition of a hidden state. These models can incorporate more than two states; we show four states above. Note than that in addition to the same speciation (λ) and extinction (μ) parameters, the MuSSE model as shown allows for anagenetic transitions between any of the states (q). In RevBayes, it is simple to control the transition matrix if your own empirical system required it. For example, if it was biologically unrealistic for an anagenetic transition to occur from state 0 to state 2, it would be possible to set the transition matrix to only allow transitions into state two from either state 1 or 3. This is discussed more in the MuHiSSE tutorial. 


In this series of tutorials, we will not be running BiSSE and MuSSE because taking the time to perform their more complex, hidden state extensions in a Bayesian framework inherently includes their more simpler versions. In other words, we can interpret results from only the A states (0A, 1A, 2A, and 3A) from the MuHiSSE if we wanted to know the results of a non-hidden state model. But by running the more complex model in RevBayes, we have the information for both and the confidence that our results will be less likely to be hampered by false positive inferences. 


## Tensorphylo
Because Bayesian SSE models are computationally intensive, especially for large trees, TensorPhylo is a program that can improve computational time. Because it might be impractical for some users to implement TensorPhylo, we have versions of this tutorial with and without Tensorphylo.

## Workflow & Data for tutorials 
The data we will be analyzing for this tutorial is the same as discussed in Zenil-Ferguson et al. (XXXXXXX): a reanalyzed dataset of restios (Polaes:Restionaceae from [Bouchenak-Khelladi and Linder 2017](https://onlinelibrary.wiley.com/doi/abs/10.1111/evo.13364)), which are a family of bamboo-like flowering plants mostly found in the Cape Floristic Region of South Africa. This study wanted to understand how two binary habtiat traits, precipitation (dry or wetlands) and elevation (costal or mountainous), contribute to faster rates of diversification.They also combined these two traits into four states for a MuSSE model.

They had two hypotheses about what contributed to higher diversification: 
 1) The older mountainous wetland habitats are more fragmented and have promoted lineage diversification because of higher rates of turnover. 
 2) The more recent coastal dryland habitats have had much more connectivity and geographic area in which to diversify.
 
The authors fitted BiSSE, HiSSE, ClaSSE, and MuSSE in the R packages hisse and diversitree, both likelihood methods, and determined the best model with the Akaike Information Criterion (AIC). When considering elevation and precipitation independently (BiSSE, HiSSE, and ClaSSE), they found no relationship with net diversification for either trait (supporting the 2-state character independent model, CID-2). When they combined the traits in a 4-state MuSSE model, the best-fit model was the one that had the same speciation but different extinction rates. Overall, Bouchenak-Khelladi and Linder (2017) concluded that no habitat category had higher net diversification rates, thus not supporting their original hypotheses. 

So how can we understand diversification dynamics at a deeper level when we implement bayesian statistics to this question? 

To start the tutorial, navigate to the appropriate folder for whether or not you want to use tensorphylo or not (i.e., "1_HiSSE" or "2_HiSSE-tensorphylo") and read the markdown file in each. The folders will have the appropriate input and output files. 



