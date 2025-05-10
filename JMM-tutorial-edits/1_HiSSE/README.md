# Bayesian Hypothesis Testing with SSE models in RevBayes 
## Tutorial overview
The purpose of this tutorial is to describe Bayesian hypothesis testing with SSE models, specifically with Hidden State Speciation and Extinction (HiSSE) and Hidden Multiple State Speciation and Extinction (MuHiSSE) models in RevBayes. 

Both HiSSE and MuHiSSE estimate differences in state-dependent diversification but differ in the number of allowed states: HiSSE requires binary states (0 and 1) and MuHiSSE allows for more than two states. HiSSE and MuHiSSE are the hidden state extensions, respectively, of the Binary Speciation and Extinction (BiSSE) and Multiple State Speciation and Extinction (MuSSE) models. It is important to mention here that all these models only allow for a trait change to occur along a branch (anagenetic trait evolution, denoted by the q transition parameters) rather than at a speciation event (cladogenetic trait evolution). If this is something you would like to implement in your own work, see the Cladogenetic Speciation and Exinction (ClaSSE, XXXX paper). 

<img width="975" alt="Screen Shot 2025-05-07 at 2 32 37 PM" src="https://github.com/user-attachments/assets/eb26d248-be98-42b3-9746-948415ea3385" />
HiSSE is a model extension of the BiSSE model by its addition of a hidden state. Both have two states, 0 and 1, and HiSSE has the additional layer of complexity by having A and B versions of the two states. Each state has its own speciation (λ, lambda) and extinction (μ, mu) parameters as well as transitions between states (q parameters). HiSSE additionally has parameters that govern the transitions between the A and B states (alpha and beta). 
<img width="1184" alt="Screen Shot 2025-05-07 at 2 58 43 PM" src="https://github.com/user-attachments/assets/256a7b8d-43be-4a4d-b547-4698d0c696a5" />


MuHiSSE is a model extension of the MuSSE model by its addition of a hidden state. These models can incorporate more than two states; we show four states above. Note than that in addition to the same speciation (λ) and extinction (μ) parameters, the MuSSE model as shown allows for anagenetic transitions between any of the states (q). In RevBayes, it is simple to control the transition matrix if your own empirical system required it. For example, if it was biologically unrealistic for an anagenetic transition to occur from state 0 to state 2, it would be possible to set the transition matrix to only allow transitions into state two from either state 1 or 3. This is discussed more in the MuHiSSE tutorial. 


In this series of tutorials, we will not be running BiSSE and MuSSE because taking the time to perform their more complex, hidden state extensions in a Bayesian framework inherently includes their more simpler versions. In other words, we can interpret results from only the A states (0A, 1A, 2A, and 3A) from the MuHiSSE if we wanted to know the results of a non-hidden state model. But by running the more complex model in RevBayes, we have the information for both and the confidence that our results will be less likely to be hampered by false positive inferences. 



Because Bayesian SSE models are computationally intensive, especially for large trees, TensorPhylo is a program that can improve computational time. Because it might be impractical for some users to implement TensorPhylo, we have versions of this tutorial with and without Tensorphylo. 
