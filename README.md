## Analysis of desialylated and sialylated spectra of Myozyme<sup>®</sup>

This repository contains R code  that was used in Myozyme<sup>®</sup> analysis, to integrate the composition of glycans in experimentally desialylated protein spectrum with _in silico_ desialylated intact glycoform annotations. 
The annotations of Myozyme<sup>®</sup> glycoforms were corrected in multiple steps:

1. The intact glycoform annotations were computationally desialylated to obtain a desialylated _in silico_ spectrum of Myozyme<sup>®</sup> 

2. The _in silico_ desialylated masses were filtered based on their correspondence with the experimentally desialylated masses

3. The _in silico_ desialylated hit scores were filtered with a cut-off (range from 0.01 to 1)

4. Relative abundances of computational and experimental Myozyme spectra were calculated and plotted

The complete analysis and plotting of the publication figures can be followed step-be-step in desialylation.qmd document which was rendered either as desialylation.html or desialylation.pdf. 