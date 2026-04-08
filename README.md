## The biodiversity of plant-frugivore interactions: types, functions, and consequences

Data and code accompanying paper in Oikos Special Issue on Frugivory and Seed Dispersal (for the FSD 2024 symposium, Ilhéus, Brazil). doi: [10.1101/2025.11.14.688423](https://doi.org/10.1101/2025.11.14.688423)

This repository includes both datasets and code for analyses. Please see the main paper and the Supplementary Material.

[Access this dataset on Dryad.](https://datadryad.org/dataset/doi:10.5061/dryad.6t1g1jxdb)
##### author:
#### name: Pedro Jordano

##### affiliation:
##### name: Estación Biológica de Doñana, CSIC and Dept. Biol. Vegetal y Ecología, Universidad de Sevilla

date: 2025-01-21

orcid: [0000-0003-2142-9116](https://orcid.org/0000-0003-2142-9116)

email: jordano@ebd.csic.es
city: Sevilla, Spain

url: [](http://pjordanolab.ebd.csic.es)

_keywords_: seed dispersal, mutualism, plant-animal interactions, frugivory

### Description of the data and file structure
***Summary***:
Pairwise plant-frugivore mutualistic interactions build up into mega-diverse networks involving dozens of interacting species, being the most generalized among free-living species. These mutualisms consist of food provisioning by plants and, their counterpart, plant propagule (seeds) movement by the animals, being crucial for the natural vegetation regeneration in many ecosystems. Yet we are far from understanding which part of this enormous interaction biodiversity is needed for their maintenance.
I overview the diversity of interaction modes involved in these mutualisms, the main components of the seed dispersal services, and their functional diversity. I examine how interaction richness covaries with partner species richness at different scales, resulting in variable patterns of species complementarities in terms of seed dispersal effects. The functionality of most generalized plant-frugivore mutualisms relies on complementarity of effects across a high diversity of partners, yet frequently depends on just a distinct subset of them, resulting in high functional redundancy. Two distinct aspects are relevant: 1) variable quantitative effects among species; 2) variable parwise-interaction outcomes, between the extremes of antagonism and mutualism.

### Sharing/Access information
Links to other publicly accessible locations of the data:
- GitHub repository: https://github.com/pedroj/oikos_2026_si_frugivory
- Data on _Euterpe edulis_ can be accessed [here](https://datadryad.org/dataset/doi:10.5061/dryad.2pm42#citations).

### Dataset contents
Folder and file structure:
- _data/data.tsv_
Biodiversity data for plant-frugivore mutualistic assemblages. 
Dataset: number of interactions and nodes in plant-frugivore interaction
networks.
Variables included:   

    ***id***: number ID. of dataset.
    ***code***: sp, species-based network; ind: individual-based network.
    ***sp***: plant species name and family.
    ***P***: number of plant species (for sp-based netwotks) or plant individuals (for individual-based networks) studied. 
    ***A***: number of animal species.
    ***S***: total number of nodes in the network (A+P).
    ***I***: number of distinct interactions in network.
    ***C***: connectance (C= i/(A+P)).
    ***biome***: t, tropical; nt: non-tropical.
    ***site***: study locality. 
    ***reference***: Bibliographic reference.

Reference for data: Quintero et al. 2025. PNAS, doi= 10.1101/2024.02.02.578595.        

-------
- _data/frug-mut_gradient.tsv_
Body mass, percent mutualistic interaction outcomes, and type of fruit consumption/dispersal for the assemblages of _Prunus mahaleb_ and _Euterpe edulis_.
From Jordano and Schupp (2000) and Galetti et al. (2013).
Variables included:
    ***species***: Frugivore species.
    ***sp_plant***: Plant species.
    ***type***: Frugivore type (seed disperser, pulp consumer, seed predator, seed disperser/pulp consumer, seed predator/pulp consumer).
    ***pct_mut***: percent fruits handled with seed ingested and dispersed away from plant.
    ***body_mass***: Body mass (g).
-------
- data/funct.csv
Functional data. _Prunus mahaleb_. _Euterpe edulis_.
Data from Jordano and Schupp 2000. Jordano et al. 2025. and Galetti et al. 2013. Additional morphological data for _Euterope edulis_ from Bovo et al. 2018. For _Prunus mahaleb_. abundance data derived from N= 43 transect censuses of total length 66.4 km and totalling 4218 birds recorded. For Euterpe edulis. 
Variables included:
***species***: Frugivore species.
***plant***: Plant species.
***type***: Frugivory type (see above).
***frug***: Percentage of diet made up of fleshy fruits (by volume).
***bmass***: Body mass (g).
***gape***: Gape width (mm) measured at the mouth commisures.
***vis.10h***: Visit rate to the tree, no. visits/10 h.
***pct_mut***: Perecentage of mutualistic interactions: based on % fruit swallowed.
***birds.km***: No. birds recorded per km census. 
***contrib_FD***: Contribution to _FD_ value. Species-specific contribution to _FD_. The product of the species' visitation rate (vis.10h) and total FD value for the plant. divided by the  total  number of visits of all frugivore species (sum(vis.10h)).
***eff_per_vis***: Seed disperser effectiveness per visit: no. seeds estimated to be removed from maternal tree per visit.
***eff_total***: vis.10h * eff_per_vis.
***prop_disp_service***: Proportion of species eff_total relative to the pooled value for the whole frugivore assemblage.

Abundance data (birds.km) refers to IPA from transect censuses in non-defaunated areas of SE Brazil Mata Atlantica.
Visitation data (***no. visits per 10h***) were obtained by direct focal watches at fruiting trees for 336 h (16 tree-days per two seasons) in _Prunus mahaleb_and 2326 h in _Euterpe edulis.

---------
- data/SDE_euterpe.txt
- data/SDE_prunus.txt
Two datasets with effectiveness components for the two studied plant species.
_Prunus mahaleb_ frugivore assemblage. Jordano & Schupp (2000) _Ecological Monographs_.
Visitation data come from 107.3 h direct watches.
***abundance***- Mean no. birds censused/km, averaged for two study years.
***visits***- Mean no. visitis recorded to fruiting trees (/10 h).
**prop_visits**- Proportion of total visits recorded (feeding records) contributed by species. Relative to the total no. records in two study years.
***eff_per_vis***- Mean no. fruits swallowed per visit (successfully dispersed seeds).
**eff_total**- Visit rate * eff_per_vis*prop fruits swallowed
***prop_disp_service***-  Proportion of total dispersal service contributed by species.

For _Euterpe edulis_ frugivore assemblage. Galetti et al. (2013) _Science_.
***species***: Frugivore species.
***behav***: Foraging and seed dissemination behaviour (defecation/regurgitation).
***mass***: Body mass (g).
***gape***: Gape width (mm).
***frugscore***: frugivory score, from 1 (sporadic fruit consumption) to 3, extremely frugivorous on _E. edulis_.
***nvis10h***: Visitation rate (/10 h).
***nfrhand***: No. fruits handled per visit.
***nfrdisp***: No. fruits ingested.
***dispprob***: Probability a seed is dispersed away from palm.
***qc***: Quantitative component of effectiveness.

---------
- data/euterpe_data.zip
  Container with source data on _Euterpe edulis_. Includes functional traits, feeding rates of animal frugivores, abundance data, and presence in different areas. Please refer to Galetti et al. 2013.
- data/prunus_data.zip
  Container with source data on _Prunus mahaleb_. Includes the data for analyses of multilayer structures in R.
- data/multiplex.zip
  Container with source data on _Prunus mahaleb_. Includes the data for analyses of multilayer structures in R: data for muxViz analyses, adjacency matrices for layers, and configuration files, as well as attributes file with individual trees characteristics.

----------
- ms_analyses_files.zip
  Container with the figures in pdf format.
----------

### Code/Software
The ***RCode*** folder (compressed file Rcode.zip) contains scripts, .Rmd, and .qmd files for data analyses in R.   
The /functions folder contains scripts with R function for multilayer analyses in muxViz and infomap.     
The /multiplex folder contains code in .Rmd and .qmd files including commented notes and thorough descriptions of the analyses. A part of the analyses, for infomap, are done with binaries included in the /multiplex/src-exe/ subfolder.  
The analysis flow is described and delineated in the files: ms_analyses.qmd, emln_analyses.qmd, infomap.Rmd, pru_muxviz.qmd, and multiplex_correlates.qmd.     
###### NOTE: Code for multilayer analysis with the muxViz package is provided in a separate subfolder (/multiplex/src-exe/). It requires specific installation of external binaries. I’ve run it on R version 4.0, since several required packages have not recent versions.

- Rcode.zip
  Container with all the code and analysis files, in R. It includes two folders, as described above:
  /functions
  /multiplex

### R Packages used
The R packages used for analyses are fully listed in the Supplementary Material file.

We used R version 4.4.1 (R Core Team 2024) and the following R packages: ade4 v. 1.7.22 (Chessel, Dufour, and Thioulouse 2004; Dray and Dufour 2007; Dray, Dufour, and Chessel 2007; Bougeard and Dray 2018; Thioulouse et al. 2018), ape v. 5.8.1 (Paradis and Schliep 2019), bayestestR v. 0.15.0 (Makowski, Ben-Shachar, and Lüdecke 2019), boot v. 1.3.31 (A. C. Davison and D. V. Hinkley 1997; Angelo Canty and B. D. Ripley 2024), correlation v. 0.8.6 (Makowski et al. 2020b, 2022), datawizard v. 1.0.0 (Patil et al. 2022), devtools v. 2.4.5 (Wickham et al. 2022), easystats v. 0.7.3 (Lüdecke et al. 2022), effectsize v. 1.0.0 (Ben-Shachar, Lüdecke, and Makowski 2020), equatiomatic v. 0.3.3 (Anderson, Heiss, and Sumners 2024), FD v. 1.0.12.3 (Laliberté and Legendre 2010; Laliberté, Legendre, and Shipley 2014), geometry v. 0.5.1 (Habel et al. 2025), ggcorrplot v. 0.1.4.1 (Kassambara 2023), gridExtra v. 2.3 (Auguie 2017), here v. 1.0.1 (Müller 2020), insight v. 1.0.1 (Lüdecke, Waggoner, and Makowski 2019), kableExtra v. 1.4.0 (Zhu 2024), knitr v. 1.49 (Xie 2014, 2015, 2024), lattice v. 0.22.6 (Sarkar 2008), lme4 v. 1.1.36 (Bates et al. 2015), Matrix v. 1.7.1 (Bates, Maechler, and Jagan 2024), modelbased v. 0.8.9 (Makowski et al. 2020a), parameters v. 0.4.1 (Lüdecke et al. 2020), performance v. 0.13.0 (Lüdecke, Ben-Shachar, et al. 2021), permute v. 0.9.7 (Simpson 2022), report v. 0.6.0 (Makowski et al. 2023), see v. 0.9.0 (Lüdecke, Patil, et al. 2021), tidyverse v. 2.0.0 (Wickham et al. 2019), usethis v. 3.1.0 (Wickham et al. 2024), vegan v. 2.6.8 (Oksanen et al. 2024).

***Package citations***
Davison, A.C. and D.V. Hinkley. 1997. Bootstrap Methods and Their Applications. Cambridge: Cambridge University Press. doi:10.1017/CBO9780511802843.
Anderson, D., A. Heiss, and J. Sumners. 2024. equatiomatic: Transform Models into LaTeX Equations. https://CRAN.R-project.org/package=equatiomatic.
Angelo C., and B.D. Ripley. 2024. boot: Bootstrap r (s-Plus) Functions. Auguie, Baptiste. 2017. gridExtra: Miscellaneous Functions for Grid Graphics. https://CRAN.R-project.org/package=gridExtra.
Bates, D., M. Mächler, B. Bolker, and S. Walker. 2015. Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software 67 (1): 1–48. https://doi.org/10.18637/jss.v067.i01.
Bates, D., M. Maechler, and M. Jagan. 2024. Matrix: Sparse and Dense Matrix Classes and Methods. https://CRAN.R-project.org/package=Matrix.
Ben-Shachar, Ma.S., D. Lüdecke, and D. Makowski. 2020. effectsize: Estimation of Effect Size Indices and Standardized Parameters. Journal of Open Source Software 5 (56): 2815. https://doi.org/10.21105/joss.02815.
Bougeard, S., and S. Dray. 2018. Supervised Multiblock Analysis in R with the ade4 Package. Journal of Statistical Software 86 (1): 1–17. https://doi.org/10.18637/jss.v086.i01.
Chessel, D., A. Dufour, and J. Thioulouse. 2004. The ade4 Package – I: One-Table Methods. R News 4 (1): 5–10. https://cran.r-project.org/doc/Rnews/.
Dray, S., and A. Dufour. 2007. The ade4 Package: Implementing the Duality Diagram for Ecologists. Journal of Statistical Software 22 (4): 1–20. https://doi.org/10.18637/jss.v022.i04.
Dray, S., A. Dufour, and D. Chessel. 2007. The ade4 Package – II: Two-Table and K-Table Methods. R News 7 (2): 47–52. https://cran.r-project.org/doc/Rnews/.
Habel, K., R. Grasman, R.B. Gramacy, P. Mozharovskyi, and D.C. Sterratt. 2025. geometry: Mesh Generation and Surface Tessellation. https://CRAN.R-project.org/package=geometry.
Kassambara, A. 2023. ggcorrplot: Visualization of a Correlation Matrix Using ggplot2. https://CRAN.R-project.org/package=ggcorrplot.
Laliberté, E., and P. Legendre. 2010. A Distance-Based Framework for Measuring Functional Diversity from Multiple Traits. Ecology 91: 299–305.
Laliberté, E., P. Legendre, and B. Shipley. 2014. FD: Measuring Functional Diversity from Multiple Traits, and Other Tools for Functional Ecology.
Lüdecke, D., M.S. Ben-Shachar, I. Patil, and Do. Makowski. 2020. Extracting, Computing and Exploring the Parameters of Statistical Models Using R. Journal of Open Source Software 5 (53): 2445. https://doi.org/10.21105/joss.02445.
Lüdecke, D., M.S. Ben-Shachar, I. Patil, P. Waggoner, and D. Makowski. 2021. performance: An R Package for Assessment, Comparison and Testing of Statistical Models. Journal of Open Source Software 6 (60): 3139. https://doi.org/10.21105/joss.03139.
Lüdecke, D., M.S. Ben-Shachar, I. Patil, B.M. Wiernik, E. Bacher, R. Thériault, and D. Makowski. 2022. easystats: Framework for Easy Statistical Modeling, Visualization, and Reporting. CRAN. https://doi.org/10.32614/CRAN.package.easystats.
Lüdecke, D., I. Patil, M.S. Ben-Shachar, B.M. Wiernik, P. Waggoner, and D. Makowski. 2021. see: An R Package for Visualizing Statistical Models. Journal of Open Source Software 6 (64): 3393. https://doi.org/10.21105/joss.03393.
Lüdecke, D., P. Waggoner, and D. Makowski. 2019. insight: A Unified Interface to Access Information from Model Objects in R. Journal of Open Source Software 4 (38): 1412. https://doi.org/10.21105/joss.01412.
Makowski, D., Ma.S. Ben-Shachar, and D. Lüdecke. 2019. bayestestR: Describing Effects and Their Uncertainty, Existence and Significance Within the Bayesian Framework. Journal of Open Source Software 4 (40): 1541. https://doi.org/10.21105/joss.01541.
Makowski, D., M.S. Ben-Shachar, I. Patil, and Da. Lüdecke. 2020a. Estimation of Model-Based Predictions, Contrasts and Means. CRAN. https://github.com/easystats/modelbased.
———. 2020b. Methods and Algorithms for Correlation Analysis in R. Journal of Open Source Software 5(51): 2306. https://doi.org/10.21105/joss.02306.
Makowski, D., D. Lüdecke, I. Patil, R. Thériault, M.S. Ben-Shachar, and B.M. Wiernik. 2023. Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption. CRAN. https://easystats.github.io/report/.
Makowski, D., B.M. Wiernik, I. Patil, D. Lüdecke, and M.S. Ben-Shachar. 2022. correlation: Methods for Correlation Analysis. https://CRAN.R-project.org/package=correlation.
Müller, K. 2020. here: A Simpler Way to Find Your Files. https://CRAN.R-project.org/package=here. Oksanen, J., G.L. Simpson, F.G. Blanchet, R. Kindt, P. Legendre, P.R. Minchin, R.B. O’Hara, et al. 2024. vegan: Community Ecology Package. https://CRAN.R-project.org/package= vegan.
Paradis, E., and K. Schliep. 2019. Ape 5.0: An Environment for Modern Phylogenetics and Evolutionary Analyses in R. Bioinformatics 35: 526–28. https://doi.org/10.1093/bioinformatics/bty633.
Patil, I., D. Makowski, M.S. Ben-Shachar, B.M. Wiernik, E. Bacher, and D. Lüdecke. 2022. datawizard: An R Package for Easy Data Preparation and Statistical Transformations. Journal of Open Source Software 7 (78): 4684. https://doi.org/10.21105/joss.04684.
R Core Team. 2024. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. https://www.R-project.org/.
Sarkar, D. 2008. Lattice: Multivariate Data Visualization with r. New York: Springer. http://lmdvr.r-forge.r-project.org.
Simpson, G.L. 2022. permute: Functions for Generating Restricted Permutations of Data. https://CRAN.R-project.org/package=permute.
Thioulouse, J., S. Dray, A. Dufour, A. Siberchicot, T. Jombart, and S. Pavoine. 2018. Multivariate Analysis of Ecological Data with ade4. Springer. https://doi.org/10.1007/978-1-4939-8850-1.
Wickham, H., M. Averick, J. Bryan, W. Chang, L. D’Agostino McGowan, R. François, G. Grolemund, et al. 2019. Welcome to the tidyverse. Journal of Open Source Software 4 (43): 1686. https://doi.org/10.21105/joss.01686.
Wickham, H., J. Bryan, M. Barrett, and A. Teucher. 2024. usethis: Automate Package and Project Setup. https://CRAN.R-project.org/package=usethis.
Wickham, H., J. Hester, W. Chang, and J. Bryan. 2022. devtools: Tools to Make Developing r Packages Easier. https://CRAN.R-project.org/package=devtools.
Xie, Y. 2014. knitr: A Comprehensive Tool for Reproducible Research in R. In Implementing Reproducible Computational Research, edited by Victoria Stodden, Friedrich Leisch, and Roger D. Peng. Chapman; Hall/CRC.
———. 2015. Dynamic Documents with R and Knitr. 2nd ed. Boca Raton, Florida: Chapman; Hall/CRC. https://yihui.org/knitr/.
———. 2024. knitr: A General-Purpose Package for Dynamic Report Generation in r. https://yihui.org/knitr/.
Zhu, H. 2024. kableExtra: Construct Complex Table with kable and Pipe Syntax. https://CRAN.R-project.org/package=kableExtra.
