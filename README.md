# BMI-CNV
A Bayesian framework for multiple genotyping platforms detection of copy number variation.
## Author
Xizhi Luo, Guoshuai Cai, Alexander C. Mclain, Christopher I. Amos, Bo Cai, Feifei Xiao.
## Description
Copy number variations (CNVs) are gains or losses of chromosomal segments and have been found to be associated with many complex human diseases. Whole-exome sequencing (WES) enables detection of CNVs with high resolution in functional protein-coding regions. However, variations in the intergenic or intragenic regions are excluded from these studies. Fortunately, in many exiting large cohorts, samples have been previously sequenced by different genotyping platforms, such as SNP array. As a result, methods for integrating multiple genotyping platforms are highly demanded for improved CNV detection. Moreover, conventional single sample-based CNV calling methods often suffer from high false discovery rate due to prominent data noise. A multi-sample strategy may reduce detection error and will be more robust to data variations. 

We developed BMI-CNV, a Bayesian Multi-sample and Integrative CNV (BMI-CNV) profiling method using existing WES and microarray data. By incorporating complementary information from multiple platforms, our method can detect CNVs with a genome-wide scale, while integrating concurrent information shared by multiple samples dramatically improve the detection performance of common CNVs. For the multi-sample integration, we identify the shared CNVs regions across samples using a Bayesian probit stick-breaking process model coupled with a Gaussian Mixture model estimation. 
## General workflow
Our method mainly focuses on CNV detection by integrating the SNP array and WES data, although it can also be naturally applied to the WES data only situation. Fig 1 shows an overview of the framework of BMI-CNV. First, WES read counts and SNP array intensities are integrated using a series of data integration procedures, including normalization, standardization, and merging. Our main algorithm consists of two main stages: Stage I uses a Bayesian PSBP method  coupled with a Gaussian mixture model-based initial data filtering to identify shared CNV regions, and Stage II as the individual CNV calling procedure. 

![workflow](Fig1.jpg)
