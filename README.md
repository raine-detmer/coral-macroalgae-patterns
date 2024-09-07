# coral-macroalgae-patterns
Code for making the figures in "A Busse balloon in the lagoon: herbivore behavior generates spatial patterns in coral reef ecosystems" by A. Raine Detmer, Scott D. Miller, Alexandra K. Dubel, Kacie Ring, Christian John, Cheryl J. Briggs, Andrew Rassweiler, and Holly V. Moeller

Model analyses:
The .m scripts labeled "FigX" are named for the figures they create and contain all code for simulating and plotting the results shown in those figures. For example, running the “Fig2.m” script will make Figure 2. The remaining .m scripts contain the functions used for simulating the models and processing model output. The "code output" folder contains the stored outputs of simulations used to make the figures (each labeled for the figure to which it corresponds). All model simulations were performed in Matlab R2023b. 

Herbivore data analyses:
The "Herbivore code" folder contains the code used for analyzing the Moorea Coral Reef LTER fish survey data (Moorea Coral Reef LTER and A. Brooks 2023). These analyses were used to calculate the relative abundances of the herbivore species included in Table 2 in Moorea. They were performed in R version 4.3.1 (R Core Team 2023). 

DOI Badge for this repository:
[![DOI](https://zenodo.org/badge/843644567.svg)](https://zenodo.org/doi/10.5281/zenodo.13730757)]

References

Moorea Coral Reef LTER and A. Brooks. 2023. MCR LTER: Coral Reef: Long-term Population and Community Dynamics: Fishes, ongoing since 2005 ver 62. Environmental Data Initiative. https://doi.org/10.6073/pasta/75644add7e7f90c568bf5045264d359a (Accessed 2024-09-06).

R Core Team (2023). R: A Language and Environment for Statistical
  Computing_. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>.

