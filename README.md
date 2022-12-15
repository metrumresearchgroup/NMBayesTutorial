# NMBayesTutorial
Supplementary code for NONMEM Bayes Tutorial

Packages have been specified in _/renv_ and should be accessible by first running in your RStudio Console window:


~~~ 

install.packages("renv")
library(renv)
renv::restore()  

~~~

For further information on using renv, please see: https://rstudio.github.io/renv/articles/renv.html#reproducibility

Directory listing:

~~~
   
   /data = simulated data used in the example runs with associated specification documents
   
   /deliv = figures and tables
   
   /model = NONMEM and mrgsolve-formatted model files (.mod)

   /script = the scripts to run the example models (.R)
   
   
