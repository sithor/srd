# srd

This is an `R` package to display categorical data as scaled rectangle diagrams, rather like Euler plots. It relies on gfortran for compilation.

The package was coded by Associate-Professor Roger Marshall, formerly from the University of Auckland, who has now retired. His email address is [rogermarshall65\@gmail.com](mailto:rogermarshall65@gmail.com){.email}.

Currently, it has been tested on linux (posit cloud) and it works, despite some Fortran warnings.

It has been tested on Apple silicon and the Fortran code throws an error.

To install, please use the following code:

`if(!require(devtools)) install.packages("devtools")`\
`devtools::install_github("sithor/srd")`

Please see examples of scaled rectangle diagrams in publications such as the following:

[Marshall, Roger J. "Scaled rectangle diagrams can be used to visualize clinical and epidemiological data." Journal of clinical epidemiology 58.10 (2005): 974-981.](https://doi.org/10.1016/j.jclinepi.2005.01.018)

[Thornley S, Sundborn G, Engelman D, Roskvist R, Pasay C, Marshall R, Long W, Dugu N, Hopoi N, Moritsuka S, McCarthy J, Morris AJ. Children's scabies survey indicates high prevalence and misdiagnosis in Auckland educational institutions. J Paediatr Child Health. 2023 Dec;59(12):1296-1303. doi: 10.1111/jpc.16512. Epub 2023 Nov 2. PMID: 37920140.](https://onlinelibrary.wiley.com/doi/epdf/10.1111/jpc.16512)
