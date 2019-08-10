PCAtools: everything Principal Components Analysis
================
Kevin Blighe
2019-08-10

Principal Components Analysis (PCA) is a very powerful technique that has wide applicability in data science, bioinformatics, and further afield. It was initially developed to analyse large volumes of data in order to tease out the differences/relationships between the logical entities being analysed. It extracts the fundamental structure of the data without the need to build any model to represent it. This 'summary' of the data is arrived at through a process of reduction that can transform the large number of variables into a lesser number that are uncorrelated (i.e. the ‘principal components'), whilst at the same time being capable of easy interpretation on the original data (Blighe, Lewis, and Lun 2018) (Blighe 2013).

*PCAtools* provides functions for data exploration via PCA, and allows the user to generate publication-ready figures. PCA is performed via *BiocSingular* (Lun 2019) - users can also identify optimal number of principal components via different metrics, such as elbow method and Horn's parallel analysis (Horn 1965) (Buja and Eyuboglu 1992), which has relevance for data reduction in single-cell RNA-seq (scRNA-seq) and high dimensional mass cytometry data.

Availability
------------

-   Release version: [PCAtools](https://www.bioconductor.org/packages/release/bioc/html/PCAtools.html)

-   Development version: [PCAtools](https://www.bioconductor.org/packages/devel/bioc/html/PCAtools.html)

-   Vignette: [PCAtools: everything Principal Components Analysis](https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html)

Bibliography
------------

Blighe, K. 2013. “Haplotype classification using copy number variation and principal components analysis.” The Open Bioinformatics Journal 7:19-24.

Blighe, K, M Lewis, and A Lun. 2018. “PCAtools: everything Principal Components Analysis.” <https://github.com/kevinblighe>.

Buja, A, and N Eyuboglu. 1992. “Remarks on Parallel Analysis.” Multivariate Behav. Res. 27, 509-40.

Horn, JL. 1965. “A rationale and test for the number of factors in factor analysis.” Psychometrika 30(2), 179-185.

Lun, A. 2019. “BiocSingular: Singular Value Decomposition for Bioconductor Packages.” R package version 1.0.0, https://github.com/LTLA/BiocSingular.
