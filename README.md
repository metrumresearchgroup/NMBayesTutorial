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
   
      /mrgsolve = mrgsolve model files
      
      /pk = model control streams
      
         1000 = Example Bayesian PopPK
         
         1000-FOCE = Example FOCE PopPK
         
         1000-7.4 = NONMEM 7.4 version with FORTRAN write statements
         
         2000 = Pediatric example
         
         2001 = Sensitivity analysis (Fix Ka)
         
         2002 = Sensitivity analysis (Variance x 10)
         
         2003 = Sensitivity analysis (Variance x 100)
         
         2004 = Sensitivity analysis (Increase location 50%)
         
         2005 = Sensitivity analysis (Decrease location 50%)

   /script = the scripts to run the example models (.R)
   
      applied-sims.R = Applied simulations for pediatric example
      
      demo-model-table.R = Making model tables
      
      dfeval.R = Inverse Wishart distributions
      
      diagnostic-templates/ = diagnostic templates for model diagnostics
      
      functions-diagnostics-npde.R = Helper script for post-hoc npde generation
      
      functions-diagnostics-rhat-ess.R = Helper script Bayesian-specific diagnostics
      
      functions-diagnostics.R = Helper script for diagnostics
      
      functions-model.R = Helper script for Bayesian model management
      
      functions-table.R = Helper script for table development
      
      model-management.R = Model management script
      
      pk-model-diagnostics-report.R = Diagnostic generating script
      
      pksens.R = Sensitivity analysis example script
      
      run-sims-npde.R = Helper script for post-hoc npde generation
      
      simout/ = Simulation output for post-hoc diagnostics
      
      tags.yaml = Model tags
   
   
