# Cohesion in collective motion

## Citation
Jadhav V, Guttal V, Masila DR. 2022 Randomness in the choice of neighbours promotes cohesion in mobile animal groups. R. Soc. Open Sci. 9: 220124, doi: [https://doi.org/10.1098/rsos.220124](https://doi.org/10.1098/rsos.220124)

This repository contains codes to generate data used for analysis in the article, “Randomness in the choice of neighbours promotes cohesion in mobile animal groups”. The codes are tested to run on Matlab version R2018b and newer versions.

You are welcome to use our codes for any of your work, with attribution -- specifically, by citing the above paper. 


## Codes for Figures 2 and 3

**Generate data of positions, velocities and networks to analyse group cohesion for given group size and parameters.**

Use Matlab code `/figures/long_sim_data_pw.m` to generate all the required data to analyse group cohesion for stochastic pairwise interactions. Positions, velocities of agents are stored in the `.mat` file named `n_pw.mat`. 

Variables are commented within the code. To reproduce the results, it is advised to run for at least T = 3500 and the number of realisations (no_it) = 15. Interaction rates, neighbourhood range (K) and all the parameters can be changed in `long_sim_data_pw.m`.

After running the simulation, use the data in `n_pw.mat` to calculate the cohesion parameter, average number of clusters and group polarisation. To do so, run the Matlab code `/figures/clus_ana_pw.m`. This file stores all the required data in the `cluster_ana_pw.mat` file. Similarly, run `/figures/alt_clus_ana.m` to calculate `expanse`, an alternate measure of group cohesion.

The folder (`figures`) also contains `avg_connections.m`, which is used to construct a network at every time window tw (called `conn_time` in `long_sim_data_pw.m`). We also use code `DBSCAN.m`, written by S. Mostapha Kalami Heris, to quantify group cohesion. 

Finally, run the code `plot_fig_2_3.m` to get figures 2A, B (and inset) of the main text.

### Simulation video

To see the collective motion of agents, run the code `simulate.m` in the folder `figures`. Make sure that you load `n_pw.mat`. Variables are defined within the code and can be changed accordingly.

## Codes for Figure 4

Generate data of positions, velocities and networks to compare group cohesion between stochastic pairwise interaction and Vicsek like averaging model for given group size and parameters.

Use Matlab code `/figures/long_sim_data_at.m` to generate all the required data to analyse group cohesion for averaging type model. Positions, velocities of agents are stored in the `.mat` file named `n_at.mat`.

Similar to codes for figure 2, variables are commented within the code, and it is advised to run for at least T = 3500 and the number of realisations (no_it) = 15 to reproduce the results.

After saving the required data, use the data in `n_at.mat` to calculate the cohesion parameter, average number of clusters and group polarisation for averaging type interaction. To get these measures run the Matlab code `/figures/clus_ana_at.m`. This file stores all the required data in the `cluster_ana_at.mat` file.    	

Finally, run the code `plot_fig_4.m` to get figure 4 of the main text.

**Note:** You also need to run `long_sim_data_pw.m` and `clus_ana_pw.m` files for the same N (group size) and K values to generate figure 4.  

### Simulation video

Similar to as explained above, to see the collective motion of agents for averaging type model, run the code `simulate.m` in the folder `figures`. Make sure that you load `n_at.mat`. Variables are defined within the code and can be changed accordingly.

## Minimal model

To simulate the minimal described in Appendix B use the codes in folder `toy_model_syn`. Use Matlab code `/toy_model_syn/long_sim_data_pw.m` to generate positions, velocities and networks to analyse group cohesion in minimal model with pairwise interaction. Positions, velocities of agents are stored in the `.mat` file named `n_pw.mat`.

Similar to above sections, variables are commented within the code, and it is advised to run for at least T = 3500 and the number of realisations (no_it) = 15 to reproduce the results.

Use data in `n_pw.mat` to calculate the cohesion parameter, average number of clusters and group polarisation for minimal model with pairwise interaction. To get these measures run the Matlab code `/toy_model_syn/clus_ana_pw.m`. This file stores all the required data in the `cluster_ana_pw.mat` file. Similarly, run `/toy_model_syn/alt_clus_ana.m` to calculate `expanse`, an alternate measure of group cohesion for minimal model with pairwise interaction.

Finally, run the code `plot_coh_para.m` to get cohesion parameter (C) as a function of K.

### Simulation video

To see the collective motion of agents for minimal model, run the code `simulate.m` in the folder `toy_model_syn`. Make sure that you load `n_pw.mat`. Variables are defined within the code and can be changed accordingly.

## Citation
You are welcome to use our codes for any of your work, with attribution -- specifically, by citing the paper: 

Jadhav V, Guttal V, Masila DR. 2022 Randomness in the choice of neighbours promotes cohesion in mobile animal groups. R. Soc. Open Sci. 9: 220124, doi: [https://doi.org/10.1098/rsos.220124](https://doi.org/10.1098/rsos.220124)

