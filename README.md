Bioalerts: a python library for the derivation of structural alerts from bioactivity and toxicity data sets
=========

Isidro Cortes-Ciriano

Journal of Cheminformatics 2016 8:13

DOI: 10.1186/s13321-016-0125-7©  Cortes-Ciriano. 2016


![alt text](https://cdn.rawgit.com/isidroc/bioalerts/master/scheme.svg "Bioalerts scheme")

Abstract
=========
Background
------
Assessing compound toxicity at early stages of the drug discovery process is a crucial task to dismiss drug candidates likely to fail in clinical trials. Screening drug candidates against structural alerts, i.e. chemical fragments associated to a toxicological response prior or after being metabolized (bioactivation), has proved a valuable approach for this task. During the last decades, diverse algorithms have been proposed for the automatic derivation of structural alerts from categorical toxicity data sets.

Results and conclusions
------
Here, the python library bioalerts is presented, which comprises functionalities for the automatic derivation of structural alerts from categorical (dichotomous), e.g. toxic/non-toxic, and continuous bioactivity data sets, e.g. KiKi or pIC50pIC50 values. The library bioalerts relies on the RDKit implementation of the circular Morgan fingerprint algorithm to compute chemical substructures, which are derived by considering radial atom neighbourhoods of increasing bond radius. In addition to the derivation of structural alerts, bioalerts provides functionalities for the calculation of unhashed (keyed) Morgan fingerprints, which can be used in predictive bioactivity modelling with the advantage of allowing for a chemically meaningful deconvolution of the chemical space. Finally, bioalerts provides functionalities for the easy visualization of the derived structural alerts.

Contents
------
Bioalerts library and documentation. The folder build contains an HyperText Markup Language (HTML) tree which documents the library bioalerts using reStructuredText (.rst) as markdown language and processed with sphinx (www.http://sphinx-doc.org/). The documentation can be browsed by opening the file index.html file in any HTML browser. The documentation of the python library RDKit can be accessed at www.rdkit.org.

Bioalerts Tutorial. The tutorial is provided in three file formats: (i) file bioalerts_tutorial.html, which can be browsed in any HTML browser, (ii) file bioalerts_tutorial, and (iii) file bioalerts_tutorial.ipynb, which can be opened with ipython notebook. This tutorial illustrates the functionalities of bioalerts on three real world data sets available at http://​github.​com/​isidroc/​bioalerts.

Details
------
Operating system(s): Platform independent

Programming language: Python

License: GNU GPL version 3

Any restrictions to use by non-academics: none
