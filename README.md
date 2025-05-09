## README
 This repository contains the code necessary to replicate the analyses presented in the manuscript:
 "Investigating tipping and its predictability in noisy environments: Evaluating the impact of temporal and species response correlation."
 
 Authors and Email: Sagar Karmakar (gitsagar111@gmail.com), Amit Samadder (math.amitsamadder18@gmail.com), and Joydev Chattopadhyay
                    (joydev@isical.ac.in)


## Study summary
 Understanding and identifying factors influencing the likelihood of sudden transitions in ecological systems is a significant area of 
 scientific research. Environmental fluctuations are particularly important as they can trigger these transitions before reaching the 
 system’s condition to a deterministic tipping point. While much focus has been made on noise-induced tipping due to uncorrelated 
 environmental noise, the impact of correlated noise on multi-species systems has been relatively overlooked. Specifically, studies 
 have neglected the correlation between species responses to environmental changes on system's susceptibility to tipping. This study 
 examines various two-species ecological models representing different interaction types in noisy environments. We reaffirm that elevated
 positive temporal autocorrelations in environmental fluctuations aggravate the chance of tipping. Conversely, our key findings suggest
 that elevated positive correlations in species responses generally delay the onset of tipping, except when the system dynamics is solely driven by positive interspecific interactions. The correlation of species responses is also critical in determining the reliability of 
 Early Warning Signals (EWSs) for predicting sudden ecological changes. Our findings highlight the importance of considering the 
 similarity between species responses to environmental variability, which significantly influence the likelihood and detectability of 
 dramatic ecological transitions.


## Directory and Script Descriptions
 Below is a description of the files and directories needed to reproduce the analyses discussed in the main text.


### Lineplot_for_different_q_maintext
 This directory contains the code to generate the plots by simulating data for each model shown in Figure 3 of the main text.
 -`Lineplot_shallow lake model.jl`: For Shallow-lake model (Figure 3 (a))
 -`Lineplot_Modified RM model.jl`: For Modified RM model (Figure 3 (b))
 -`Lineplot_General two species mutualism model.jl`: For General two species mutualism model (Figure 3 (c))

### Contour_Plot_maintext
 This directory contains the code to generate the contour plots by simulating data for each model shown in Figure 4 of the main text.
 -`Contour plot_shallow_lake model.jl`: For Shallow-lake model (Figure 4 (a))
 -`Contour plot Modified RM model.jl`: For Modified RM model (Figure 4 (b))
 -`Contour plot_General two species mutualism model.jl`: For General two species mutualism model (Figure 4 (c))


### Early_warning
 (1)-`EWS_pre-transition_data.jl`:generates the pre-transition time series data (for all models) analyzed to compute the Kendall tau values, which are used to evaluate the performance of early warning signals.
 (2)-`Kendall's tau generation.r` :This R script generates the Kendall tau values for each indicator using the pre-trasition time series
 data and the script `Box_plot.r` generate corresponding boxplots shown in Figure 5 (the shallow-lake model is given as an example, but the same code and procedure is also applicable for other models).

### EWS_CSV_files
 This directory contains all CSV files with pre-transition time series data for each model, including those from the main text and Appendix C of the SI. The data is generated using the Julia script `EWS_pre-transition_data.jl`.  

 Using this dataset, one can directly recreate Figure 5 from the main text and the supplementary figures (for EWSs) in the SI using the R scripts `Kendall's_tau_generation.r` and `Box_plot.r`.  

 Each model is stored in a separate file, with each file containing 10 CSV datasets of 150 pre-transition time series arranged in columns—five datasets for each of the two species, corresponding to different values of 𝑞.

### Contour_Plot_SI
 This directory contains the code to generate the contour plots for the additional models described in Appendix C of the SI.
 
### Lineplot_for_different_q_SI
 This directory contains the code to generate the line plots for the additional models described in Appendix C of the SI.

### Fluctuation_size
 - `Fluctuation_size_data_generation.jl`: Generates the simulation data of fluctuation size.
 - `Fluctuation_size_plot.jl`: Generates the fluctuation size plots shown in Appendix F.


## Software Requirements
 - All simulations are performed using `Julia version 1.11.0`.
 - The generation of Kendall tau values and the corresponding box plots are done using `R version 3.6.3`.
 (on linux X86_64 operating system)
 



## Dependencies

 the following package versions were used:

 - - Julia - -

 Distributed v0.11.0 
 DataFrames v1.7.0
 VectorizedStatistics v0.5.10
 Plots v0.40.8
 StatsPlots v0.15.5
 Statistics v1.11.1
 LaTeXStrings v0.4.32
 CSV v0.10.15

 - - R - -

 ggplot2 v3.5.1
 reshape2 v1.4.4
 vroom v1.6.5
 foreach v1.5.2
 doParallel v1.0.17

