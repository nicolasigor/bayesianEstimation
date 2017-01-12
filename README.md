# bayesianEstimation
A practical study on bayesian estimation using Kalman Filter (Linear and Extended) and Particle Filter (SIR)



This is the result of a workshop I did in the context of bayesian estimation. The objective was to implement some of the algorithms I learned in practical examples. The codes are written in Matlab, and I used them to extract the results shown in the final presentation, which corresponds to the file "Bayesian Estimation Slides.pdf".

kalmanFilterExample.m is, as the name suggests, a linear kalman filter example that allowed me to study the effect of the initial conditions.

basicNonlinearExample.m is an extended kalman filter example to study similar effects.

highlyNonlinearExample.m implements a typical example of a highly nonlinear and non stationary system, and shows a filtering stage using both the extended kalman filter and the particle filter (SIR specifically).

batteryExample.m and pola2015.m are implementation of the extended kalman filter and the particle filter using discharge curves of Li-Ion batteries in order to estimate the actual state of charge. The model used, and the data, are the ones discussed in the paper of Pola, 2015, "Particle-Filtering-Based Discharge Time Prognosis for Lithium-Ion Batteries With a Statistical Characterization of Use Profiles", which I also uploaded. The first one, batteryExample.m, attemps to filter the data simulated from a certain model with known parameters, using as input one of the dataset extracted from the paper. The second, pola2015.m, tries to replicate the first part of the paper in which the parameters of the model are unknown and have to be estimated using the methodology presented on a training data set, and then the filter is executed on a validation data set.

The data extracted from the paper of Pola, 2015, correspond to the .csv files. The data were extracted just using a software to create approximate data based on an image.
