# Cohesion in collective motion

This repository contains codes to generate data used for analysis in the article, “Randomness in the choice of neighbours promotes cohesion in mobile animal groups”.

## Codes for Figures 2 and 3

**Generate data of positions, velocities and networks to analyse group cohesion for given group size and parameters.**

Use Matlab code `/figures/long_sim_data_pw.m` to generate all the required data to analyse group cohesion for stochastic pairwise interactions. Positions, velocities of agents are stored in the `.mat` file named `n_pw.mat`. 

Variables are commented within the code. To reproduce the results, it is advised to run for at least T = 1500 and the number of realisations (no_it) = 15. Interaction rates, neighbourhood range (K) and all the parameters can be changed in `long_sim_data_pw.m`.

After running the simulation, use the data in `n_pw.mat` to calculate the cohesion parameter, average number of clusters and group polarisation. To do so, run the Matlab code `/figures/clus_ana_pw.m`. This file stores all the required data in the `cluster_ana_pw.mat` file. 

The folder (`figures`) also contains `avg_connections.m`, which is used to construct a network at every time window tw (called `conn_time` in `long_sim_data_pw.m`). We also use code `DBSCAN.m`, written by S. Mostapha Kalami Heris, to quantify group cohesion. 

Finally, run the code `plot_fig_2_3.m` to get figures 2A, B (and inset) of the main text.

### Simulation video

To see the collective motion of agents, run the code `simulate.m` in the folder `figures`. Make sure that you load `n_pw.mat` in `line 8`. Variables are defined within the code and can be changed accordingly.

## Codes for Figure 4

Generate data of positions, velocities and networks to compare group cohesion between stochastic pairwise interaction and Vicsek like averaging model for given group size and parameters.

Use Matlab code `/figures/long_sim_data_at.m` to generate all the required data to analyse group cohesion for averaging type model. Positions, velocities of agents are stored in the `.mat` file named `n_at.mat`.

Similar to codes for figure 2, variables are commented within the code, and it is advised to run for at least T = 1500 and the number of realisations (no_it) = 15 to reproduce the results.

After saving the required data, use the data in `n_at.mat` to calculate the cohesion parameter, average number of clusters and group polarisation for averaging type interaction. To get these measures run the Matlab code `/figures/clus_ana_at.m`. This file stores all the required data in the `cluster_ana_at.mat` file. 	

Finally, run the code `plot_fig_4.m` to get figure 4 of the main text.

**Note:** You also need to run `long_sim_data_pw.m` and `clus_ana_pw.m` files for the same N (group size) and K values to generate figure 4.  

### Simulation video

Similar to as explained above, to see the collective motion of agents for averaging type model, run the code `simulate.m` in the folder `figures`. Make sure that you load `n_at.mat` in `line 8`. Variables are defined within the code and can be changed accordingly.

