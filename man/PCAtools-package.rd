\name{PCAtools-package}

\title{PCAtools: everything Principal Component Analysis}

\description{Principal Component Analysis (PCA) is a very powerful technique that has wide applicability in data science, bioinformatics, and further afield. It was initially developed to analyse large volumes of data in order to tease out the differences/relationships between the logical entities being analysed. It extracts the fundamental structure of the data without the need to build any model to represent it. This 'summary' of the data is arrived at through a process of reduction that can transform the large number of variables into a lesser number that are uncorrelated (i.e. the ‘principal components'), whilst at the same time being capable of easy interpretation on the original data [@PCAtools] [@BligheK]. *PCAtools* provides functions for data exploration via PCA, and allows the user to generate publication-ready figures. PCA is performed via *BiocSingular* [@Lun] - users can also identify optimal number of principal component via different metrics, such as elbow method and Horn's parallel analysis [@Horn] [@Buja], which has relevance for data reduction in single-cell RNA-seq (scRNA-seq) and high dimensional mass cytometry data.}

\docType{package}
