Robust-Skew-t
=============

Repository to store research code for the robust skew-t fitting project.  The code will (eventually) have the following structure:

- Distribution Functions:
  - Normal
  - Skew Normal
  - Student's t
  - Skew t
  - Multivariate Normal
  - Multivariate Skew Normal
  - Multivariate t
  - Multivariate Skew t
For each distribution, we have one script which defines all the required functions: deviance, partial derivative of the deviance

- Parameter Estimation
  - Maximum Likelihood
  - Method of Moments
  - Maximum Constrained Likelihood (MCLE)
  - Trimmed Likelihood (TLE)
Each of these functions should take a distribution object defining the distribution (deviance, derivative of deviance, etc.) and a dataset.  Additionally, they may take other parameters (such as the level of trimming or initial estimate) and will return estimates of the parameters.

- Additional functions
  - Dampening Functions for MCLE
