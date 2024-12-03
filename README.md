# MMFPCA-draft

## Contents of this folder

- mmfpca_tutorial.R: A step-by-step implementation of the proposed multilevel multivariate FPCA and the associated procedures including simulated data generation and model fitting described in "Joint Modeling of Evoked and Induced Event-Related Spectral Perturbations."
  
- eigenfunction_construction.R: Functions constructing the multilevel multi- and uni-varaite two-dimensional eigenfunctions described in the Section 3 of "Joint Modeling of Evoked and Induced Event-Related Spectral Perturbations."
  
- simulation_data_generation.R: Functions simulating multilevel multivariate two-dimensional functional data under the simulation designs stated in Section 3 of "Joint Modeling of Evoked and Induced Event-Related Spectral Perturbations."
  
- mmfpca.R: Function of the algorithm fitting multilevel multivariate FPCA proposed in Section 2 of "Joint Modeling of Evoked and Induced Event-Related Spectral Perturbations."

## Introduction

The contents of this folder allow for the implementation of the estimation algorithm for fitting multilevel multivariate FPCA for high dimensional functional data proposed in "Joint Modeling of Evoked and Induced Event-Related Spectral Perturbations". Users can simulate a two-level bivaraite two-dimensional functional sample dataset and apply the proposed estimation algorithm to fit the multilevel multivariate FPCA model. Detailed instructions on how to perform the aforementioned procedures are included in mmfpca_tutorial.R.

## Requirements

The included R programs require R 4.3.1 and the packages listed in mmfpca_tutorial.R.

## Installation

Load the R program files into the global environment and install the required packages using commands in mmfpca_tutorial.R.
