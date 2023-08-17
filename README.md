# local-dependence-toolbox
This Matlab toolbox helps readers to reproduce sections 3,4,5 in the paper "Generalized Local Kendall’s Tau: A Nonlinear Local Dependence Framework".

# Description
1. Folder
   "Data": data.
   "Estimation of mixture copulas" : codes for estimating Archimedean mixture copula models and Patton's SJC copula.
   "Simulation study" : codes for reproducing the simulation results. 

2. File description
   main.m : estimating local Kendall's tau between Microsoft and IBM, and its bootstrap standard error and 95% confidence interval.
   main_diagonal_local_tau.m : Type II empirical and theoretical local Kendall's tau along the main diagonal.
   minor_diagonal_local_tau.m : Type II empirical and theoretical local Kendall's tau along the minor diagonal.
   type_I_local_tau_surface.m : Type I local Kendall's tau surface.
   type_II_local_tau_surface.m : Type II local Kendall's tau surface.
   fun_ldcopula_general.m : function for calculating Type I local Kendall's tau of the copula model.
   fun_ldcopula_type_II : function for calculating Type II local Kendall's tau of the copula model.
   fun_ldcopulasurf_general.m : function for drawing Type I local Kendall's tau surface of a copula.
   fun_ldcopulasurf_type_II.m : function for drawing Type II local Kendall's tau surface of a copula.
   fun_sampleld_general.m : function for calculating Type I local Kendall's tau in certain local region for sample data.
   fun_sampleld_type_II.m : function for calculating Type II local Kendall's tau in certain local region for sample data.
   fun_sampleldsurf_general.m : function for drawing Type I local Kendall's tau surface for sample data.
   fun_sampleldsurf_type_II.m : function for drawing Type II local Kendall's tau surface for sample data.
   fun_u_statistic_based_ld_general.m : function for calculating Type I local Kendall's tau based on U-statistic.
   fun_u_statistic_based_ld_type_II.m : function for calculating Type II local Kendall's tau based on U-statistic.

# Author Information
   Author: Zaixin Huang
   Email: eric.huangzaixin@gmail.com
   This version: 2023.01.01
   Notes: This toolbox is only free for scholars and students who use it for research purposes and cannot be used for commercial purposes.



