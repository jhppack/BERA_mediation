# BERA_mediation
Mediation analysis in Bayesian extended redundancy analysis with mixed outcome variables

This project aims to perform extended redundancy analysis from a Bayesian perspective by considering mediators in the relationship between multiple correlated predictors and continuous and ordinal outcomes. This project is referred to as BERA-mediation. It provides an R function (BERA.med) to conduct this analysis and an R function (BERA.med.summary) to summarize and organize the analysis output. Additionally, the project includes an example R code that applies these functions to Simulation Setting 1 of a currently under-review paper.

In the current version of the R function, there are some important considerations:
- The number of ordinal outcomes must be at least two.
- The product of the component weight matrix W and the regression coefficient matrix A2 (or A3) is inherently sign ambiguous, so caution is required when interpreting the results.
