# coral-macroalgae-patterns
Code for making the figures in "A Busse balloon in the lagoon: herbivore behavior generates spatial patterns in coral reef ecosystems" by A. Raine Detmer, Scott D. Miller, Alexandra K. Dubel, Kacie Ring, Christian John, Cheryl J. Briggs, Andrew Rassweiler, and Holly V. Moeller

Model analyses:
The .m scripts labeled "FigX" or "FigSX" are named for the figures they create and contain all code for simulating and plotting the results shown in those figures. For example, running the “Fig1FigS1.m” script will make Figure 1 and Supplemental Figure 1. "Rev" scripts make figures used only for responses to reviewers. The "code output" folder contains the stored outputs of simulations used to make the figures (each labeled for the figure to which it corresponds; for example, "Fig2.mat" contains the outputs from the simulations run in "Fig2.m").

The "Briggs___.m" scripts contain the functions for PDE models that each use the Briggs et al. 2018 model to describe local benthic dynamics:

• "BriggsHrPDEextH.m": contains the function for running the PDE model analyzed in the main text (i.e., the version with dynamic herbivores and external herbivore recruitment)

• "Briggs2HrPDE.m": contains the function for running the version of the PDE model with two herbivore populations

• "BriggsHPDE.m": contains the function for running the version of the PDE model with no herbivore dynamics

• "BriggsHrPDEHIC.m": contains the function for running the version of the PDE model with step-wise (rather than uniform) initial herbivore distributions

• "BriggsHrPDEStoch.m": contains the function for running the version of the PDE model with stochastic external recruitment

The "altPDE.m" and "MumbyHPDE.m" scripts contain the functions for running PDE models based on two alternative models of local benthic dynamics. The "altPDE.m" function uses a model based on the model published by van de Leemput et al. (2016), while the "MumbyHPDE.m" function uses the original Mumby et al. (2007) model.   

The functions for processing model output are listed below:

• "peakfun.m": calculates metrics describing the spatial patterns (e.g., peak widths, heights, wavelengths)

• "peakfun2.m": calculates the maximum number of peaks in a given spatial pattern

• "stepfun.m": function for generating step-wise initial conditions

• "tpfun.m": function for calculating the lower tipping point (boundary of the region of bistability) for the Briggs et al. (2018) model as a function of model parameters

• "tpfunExtH.m": function for calculating the lower tipping point (boundary of the region of bistability) for the Briggs et al. (2018) model as a function of external herbivore recruitment rate

All model simulations were performed in Matlab R2023b. 

Herbivore data analyses:
The "Herbivore code" folder contains the code used for analyzing the Moorea Coral Reef LTER fish survey data (Moorea Coral Reef LTER and A. Brooks 2023). These analyses were used to calculate the relative abundances of the herbivore species in Moorea included in Table 2 of Detmer et al. (2025). They were performed in R version 4.3.1 (R Core Team 2023). 

DOI Badge for this repository:
[![DOI](https://zenodo.org/badge/843644567.svg)](https://zenodo.org/doi/10.5281/zenodo.13730757)

References

Briggs, C.J., Adam, T.C., Holbrook, S.J. & Schmitt, R.J. (2018). Macroalgae size refuge from 
herbivory promotes alternative stable states on coral reefs. PLoS ONE, 13, e0202273.

Moorea Coral Reef LTER and A. Brooks. 2023. MCR LTER: Coral Reef: Long-term Population and Community Dynamics: Fishes, ongoing since 2005 ver 62. Environmental Data Initiative. https://doi.org/10.6073/pasta/75644add7e7f90c568bf5045264d359a (Accessed 2024-09-06).

Mumby, P.J., Hastings, A. & Edwards, H.J. (2007). Thresholds and the resilience of Caribbean 
coral reefs. Nature, 450, 98–101.

R Core Team (2023). R: A Language and Environment for Statistical
  Computing_. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>.

Van De Leemput, I.A., Hughes, T.P., Van Nes, E.H. & Scheffer, M. (2016). Multiple feedbacks 
and the prevalence of alternate stable states on coral reefs. Coral Reefs, 35, 857–865.

