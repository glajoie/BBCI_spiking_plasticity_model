Guillaume Lajoie September 2016 (http://faculty.washington.edu/glajoie/wordpress/)

DISCLAIMER: This software was design for specialized scientific computation only and is not meant as a deployable product. Some parameter settings may lead to errors.

#####################################################
More details about the network model can be found in:

-Correlation-based model of artificially induced plasticity in motor cortex by a bilateralBrain-Machine Interface, Guillaume Lajoie, Nedialko Krouchev, John F. Kalaska, Adrienne Fairhall, Eberhard E. Fetz, under review, (2016), preprint: https://arxiv.org/abs/1609.00762
#####################################################

OVERVIEW: MATLAB code to perform spiking simulations of spiking network

-There are two main types of codes: (1) Full spiking simulations (2) Analytic model simulations

-There are functions that both codes call that reside in "/functions/". This folder is added to the path in each script header.

(1) Full Spiking simulations

"Spiking_simulator.m" 
-script that performs the full network simulations. 
-A switch control wether or not BBCI is turned on.
-spits out and saves data containing synaptic trajectories and some spikes.

"Spiking_plotter.m"
-script that uses data from script above and plots a bunch of figures.

NOTE: included with the codes are data files from previously ran simulations that can be plotted using  the plotting script.

(2) Analytic estimates codes

NOTE: both codes below use manufactured external cross correlations C_hat that are gaussian-shaped, as in the paper.

"Analytic_equilibria_finder.m"
-Script that takes the same parameter as a full spiking simulation, then computes the mean synaptic equilibria with and without BBCI. 
-Script also produces plots of before and after.

"Analytic_para_swipe.m"
-Script that uses the analytic model and computes with and without BBCI equilibria for a range of parameters.
-Parameters to be explored include: C_hat GaussianSD, stimulation delay, axonal delay, dendritic delay
