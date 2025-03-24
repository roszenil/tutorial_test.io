---
layout: workshop
type: workshop
title: Bayesian Phylogenetics in RevBayes
location: Workshop on Molecular Evolution at MBL, Woods Hole, MA
startdate: 2018-07-20
enddate: 2018-07-29
registration: https://molevol.mbl.edu

description: 
- <p>Instruction in Bayesian phylogenetic inference and divergence-time estimation will be taught at <a href="https://molevol.mbl.edu/index.php/Main_Page">the Workshop on Molecular Evolution</a> at the <a href="http://www.mbl.edu/">Marine Biological Laboratory</a> (MBL). This course was founded in 1988 and is the longest-running workshop serving the field of evolutionary biology. Students work closely with internationally-recognized scientists, receiving (i) high-level instruction in the principles of molecular evolution and evolutionary genomics, (ii) advanced training in statistical methods best suited to modern datasets, and (iii) hands-on experience with the latest software tools (often from the authors of the programs they are using). The material is delivered via lectures, discussions, and bioinformatic exercises motivated by contemporary topics in molecular evolution. A hallmark of this workshop is the direct interaction between students and field-leading scientists. The workshop serves graduate students, postdocs, and established faculty from around the world seeking to apply the principles of molecular evolution to questions of anthropology, conservation genetics, development, behavior, physiology, and ecology. The workshop also welcomes participants from federal agencies and science journalists. A priority of this workshop is to foster an environment where students can learn from each other as well from the course faculty.</p><p>For the full workshop content, list of faculty, and schedule, please see the <a href="https://molevol.mbl.edu/index.php/Schedule">main course website</a>.</p><p>Instructions for working with the RevBayes tutorials on MBL computing resources are provided in the section on <a href="woodshole2018#using-revbayes-on-the-mbl-cluster">Using RevBayes on the MBL Cluster</a> below.</p>

instructors:
- <a href="http://mlandis.github.io/">Michael Landis</a>
- <a href="http://phyloworks.org/">Tracy Heath</a>
- <a href="https://ib.berkeley.edu/people/faculty/huelsenbeckj">John Huelsenbeck</a>
- <a href="https://phylogeny.uconn.edu/">Paul Lewis</a>

software:
- <a href="https://revbayes.github.io/revbayes-site/software">RevBayes v1.0.8</a> 
- <a href="http://tree.bio.ed.ac.uk/software/tracer/">Tracer</a>
- <a href="http://tree.bio.ed.ac.uk/software/figtree/">FigTree</a>
- A simple text editor, such as Sublime Text, NotePad++, TextWrangler, BBEdit, vim, or emacs 

schedule:
  - starttime: 2018-07-22T09:00
    endtime: 2018-07-22T12:00
    topic: "Lecture: An Introduction to Bayesian Phylogenetics"
    material: <a href="https://molevol.mbl.edu/index.php/Paul_Lewis">Lecture slides and other materials</a>
    instructors: Paul Lewis
  
  - starttime: 2018-07-22T12:00
    endtime: 2018-07-22T14:00
    topic: Lunch
  
  - starttime: 2018-07-22T14:00
    endtime: 2018-07-22T17:00
    topic: "Lecture: Likelihood-based Phylogenetic Inference"
    material: <a href="https://molevol.mbl.edu/index.php/John_Huelsenbeck">Lecture slides and other materials</a> 
    instructors: John Huelsenbeck

  - starttime: 2018-07-22T17:00
    endtime: 2018-07-22T19:00
    topic: Dinner

  - starttime: 2018-07-22T19:00
    endtime: 2018-07-22T20:30
    topic: "Lecture: Phylogenetic Inference and Graphical Models"
    material: <a href="https://molevol.mbl.edu/images/c/ce/Phylo_pgm_wh_2018_mlandis.pdf">Lecture slides</a>
    instructors: Michael Landis

  - starttime: 2018-07-22T20:30
    endtime: 2018-07-22T20:45
    topic: Break

  - starttime: 2018-07-22T20:45
    endtime: 2018-07-22T22:00
    topic: "Lecture: Bayesian Divergence-Time Estimation"
    material: <a href="https://figshare.com/articles/Bayesian_Divergence-Time_Estimation_Lecture/6849005">Lecture slides</a>
    instructors: Tracy Heath

  - starttime: 2018-07-23T09:00
    endtime: 2018-07-23T22:00
    topic: <a href="https://molevol.mbl.edu/index.php/Schedule">See course schedule</a>

  - starttime: 2018-07-24T09:00
    endtime: 2018-07-24T22:00
    topic: <a href="https://molevol.mbl.edu/index.php/Schedule">See course schedule</a>

  - starttime: 2018-07-25T09:00
    endtime: 2018-07-25T22:00
    topic: Free Day

  - starttime: 2018-07-26T09:00
    endtime: 2018-07-26T22:00
    topic: <a href="https://molevol.mbl.edu/index.php/Schedule">See course schedule</a>

  - starttime: 2018-07-27T09:00
    endtime: 2018-07-27T22:00
    topic: <a href="https://molevol.mbl.edu/index.php/Schedule">See course schedule</a>

  - starttime: 2018-07-28T09:00
    endtime: 2018-07-28T9:30
    topic: "Lab: Introduction to RevBayes and the Rev Language"
    material: <a href="https://molevol.mbl.edu/images/3/32/Lab_phylo_pgm_mlandis_WH2018.pdf">Lab slides</a>
    instructors: Michael Landis

  - starttime: 2018-07-28T09:30
    endtime: 2018-07-28T10:00
    topic: "Lab: Introduction to RevBayes and the Rev Language"
    tutorials: intro
    instructors: Michael Landis

  - starttime: 2018-07-28T10:00
    endtime: 2018-07-28T12:00
    topic: "Lab: Phylogenetic Inference in RevBayes"
    tutorials: ctmc
    instructors: Michael Landis
  
  - starttime: 2018-07-28T12:00
    endtime: 2018-07-28T14:00
    topic: Lunch
  
  - starttime: 2018-07-28T14:00
    endtime: 2018-07-28T14:20
    topic: "Lecture: Bayesian Divergence-Time Estimation, continued"
    material: <a href="https://figshare.com/articles/Bayesian_Divergence-Time_Estimation_Lecture/6849005">Lecture slides</a>
    instructors: Tracy Heath

  - starttime: 2018-07-28T14:20
    endtime: 2018-07-28T17:00
    topic: "Lab: Divergence-Time Estimation in RevBayes"
    tutorials: fbd/fbd_specimen
    instructors: Tracy Heath
---

## Using RevBayes on the MBL Cluster 

### Executing RevBayes

For the tutorials used in this workshop, you may want to use the version of RevBayes on the MBL cluster. 

Once you log in, you can use RevBayes from any directory on the cluster by just typing `rb`.

```
rb
```
{:.bash}


```
Build from master (94149b) on Mon Jul 23 15:55:28 EDT 2018

Visit the website www.RevBayes.com for more information about RevBayes.

RevBayes is free software released under the GPL license, version 3. Type 'license()' for details.

To quit RevBayes type 'quit()' or 'q()'.


>
```
{:.Rev-output}

### Downloading Data Files

Some of the tutorials you will us have data files associated with them. You will find hyperlinks to all of the required files on the tutorial page. These should be listed in the "Data files and scripts" box at the top of the tutorial page below the table of contents.

For this workshop, we recommend that you create a new folder for each tutorial that you work through. 
Follow the directions in the tutorial for creating your directory structure. 
Then download the individual files using `wget`.

For example, in the tutorial on divergence times, you will need to download the sequence data file:

```
wget https://revbayes.github.io/tutorials/fbd/data/bears_cytb.nex
```
{:.bash}

### Writing Rev Scripts

Some tutorials require that you write several modular Rev scripts to define your model and MCMC analysis. 
We reccomend that you compose these files on your own computer or using the text editor on the MBL cluster.
If you create these files on your personal machine, you will have to upload them to the cluster using `scp` or cyberduck.





