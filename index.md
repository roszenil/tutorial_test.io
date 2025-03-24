---
---

<div class="row" style="display: flex; align-items: center;">
<div class="col-sm-2 col-xs-4">
  <img src="{{ site.baseurl }}{% link assets/img/gm-tree.svg %}" alt="" width="100%" />
</div>
<div class="col-sm-10 col-xs-8">
  <h1 style="font-family:'Raleway'; font-size: 2.5em"><font color="#2D28B5"><b>RevBayes</b></font></h1>
  <h3>Bayesian phylogenetic inference using probabilistic graphical models and an interpreted language</h3>
</div>
</div>


### About

RevBayes provides an interactive environment for statistical computation in phylogenetics. It is primarily intended for modeling, simulation, and Bayesian inference in evolutionary biology, particularly phylogenetics. However, the environment is quite general and can be useful for many complex modeling tasks.

RevBayes uses its own language, Rev, which is a probabilistic programming language like [JAGS](http://mcmc-jags.sourceforge.net/), [STAN](http://mc-stan.org/), [Edward](http://edwardlib.org/), [PyMC3](http://docs.pymc.io/intro.html), and related software. However, phylogenetic models require inference machinery and distributions that are unavailable in these other tools.

The Rev language is similar to the language used in R. Like the R language, Rev is designed to support interactive analysis. It supports both functional and procedural programming models, and makes a clear distinction between the two. Rev is also more strongly typed than R.

RevBayes is a collaboratively [developed]({% page_url developer %}) software project. 
