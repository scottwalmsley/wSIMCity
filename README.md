#! As we refine the usability of our code, please feel free to contribute by submitting issues! Your insight is always welcome !#
# wSIMCity:</br></br>Searching for DNA-adducts in DIA wide-SIM/MS<sup>2</sup> data

</br>

## Table of Contents  


[Background](#background)

[Overview of software](#overview)</br>

[Scoring](#scoring)

[Getting started](#getstart)

[Installation](#installation)</br>

[Dependencies](#dependencies)</br>

[Before you start](#beforestart)</br>

[Usage](#usage)</br>

</br></br>

## Background




Welcome to wSIMCity, software for mining Wide Sim mass spectrometry data.  wSIMCity was developed by Scott Walmsley, PhD, of the Masonic Cancer Center and the Institute for Health Informatics at the Univerity of Minnesota - Twin Cities for the laboratory of Prof. Rob Turesky, PhD.  Development of wSIMCity was brought about by necessity to overcome data structure issues in wide-SIM data independent acquisition (DIA) data produced in DNA-adduct mass spectrometry (MS) experiments, and to facilitate automated detection of DNA-adducts.  As such we developed an R-package to process raw wide-SIM DIA data and to mine for DNA-adducts in that data.  

Wide-SIM data independent acquisition (DIA) methods were developed to increase sensitivity for detection of DNA-adducts in digests of DNA and using mass spectrometry. More specfically, the automatic gain control (AGC) of orbitrap MS is used together with wide selected ion monitoring (wide-SIM) *m/z* ranges to boost sensitivity for the ions of interest. The wide-SIM MS<sup>1</sup> data is collected and then the subsequent scan undergoes a low energy collision induced dissociation (CID) at ~25 millielectron volts (meV).  This MS<sup>2</sup> scan leverages the known neutral loss at low energy CID of a dexoyribose and serves as the tell tale sign of the presence of a DNA-adduct.   Turesky and group pioneered this altered DIA methodology and it likely will serve as a strategy for the untargeted discovery of novel DNA-adducts. The general mechanism is shown below using the molecule dG-C8-PhIP, an important carcinogen of prostate cells caused by ingestion of overcooked meats:</br>


<p align="center">
<img src="images/mechanism.png" title = "mechanism">
</p>

The blue portion of the molecule is the mutagen or carcinogen, shown bound to the nucleoside (black and red).   MS<sup>2</sup> CID (25meV) causes the neutral loss of a deoxyribose (dR, red) leading to the formation of Guanine-C8-PhIP.  Generally speaking, it has been shown that this mechanism is universal across as wide spectrum of DNA adducts.   The ions at the MS<sup>1</sup> level are known as presursor ions (denoted: [M+H]<sup>+</sup>), and the resulting ion after the neutral loss is called the 'aglycone' (denoted as: [B+H<sub>2</sub>]<sup>+</sup>) because of the neutral loss of a deoxyribose.  Note the gain of a mobile proton on the aglycone molecule, hence the '+H<sub>2</sub>' notation.

Targeted extraction of molecules such as dG-C8-PhIP from this type of data is fairly straightforward.   However, a major goal is to facilitate untargeted detection of DNA-adducts using this DIA based method. Below is a figure showing the slight differences between scanning methods common place in proteomics and metabolommics (A), and methods developed for DNA adductomics (B and C). Included in panel B is a DDA method for DNA adductomics  called constant neutral loss (CNL) screening.  Panel C is the method for which this software was originally written.  It is the very nature of the data structure produced in C that lead to incompatibilities with current DDA /DIA software developed for proteomics and metabolomics.

</br></br>

<p align="center">
<img src="images/DIA_scanning.png" title = "scanning" width ="600">
</p>

</br></br></br>

<a name="overview"/>

## Overview of wSIMCity

wSIMCity seeks to discover the 'landscape' or map of DNA-adducts in a DNA sample prepared for and analyzed using wide-SIM-MS.  The overall workflow is shown in the next figure:

</br></br>

<p align="center">
<img src="images/total_model.png" title = "scanning" width ="600">
</p>


</br></br></br>


Software developed for the DIA methods (and previously data dynamic acquisition (DDA) methodologies) have worked well, but are incompatible with our current wide-SIM scanning technique.   Therefor we have implemented a simple workaround using the R package ```mzR```, incorporated the MSDIAL feature finding algortihm, and devised a custom global modelling and search strategy to detect the adducts.   Extracted ion chromatograms (EIC) for candidate molecules are produced.  Further we devised a scheme to best levereage a global scoring method together with an individual scoring metric to help the researcher identify and pursue candidate DNA-adduct molecules for identification.   The output of wSIMCity is list of masses and retention times to be used for targetted identification using high energy dissociation (HCD) MS. 

## Scoring<a name="scoring"/>

The best aspect of wSIMCIty is it produces a series of scores to help the researcher rummage through the list of candidate DNA-adducts.  

The score is broken down into several components. The first components are borrowed from MSDIAL's scoring method for matching detected compounds with those listed in a mass and retention database.  Instead, we use this scoreing method to match precursor - aglycone molecules found in DNA-adductomics data. 

A major difference between our scoring method and MSDIAL's is that the underlying assumptions about the distributive properties of the measured values in our data are a little different:  we use ppm mass errors to differentiate between real [M+H]<sup>+</sup> and [B+H<sub>2</sub>]<sup>+</sup> pairs of features and false hits. Specifically, for every [M+H]<sup>+</sup> ion, we compute a theroetical [M+H-dR]<sup>+</sup> *m/z*, and then look for [B+H<sub>2</sub>]<sup>+</sup> ions with that mass and retention time. We look for every ion with that *m/z* value and then measure how far off in ppm the experimental ion is. This measurement is called 'Mass err (ppm)' in the below plot.  

These measured data points follow a Laplace distribution (well, it's more like a Cauchy distribution, but that is near imposible to deal with in terms of the math involved). We also, like MSDIAL, incorporate retention time as a metric.  Like MSDIAL, our assumption here follows a Guassian distribution. [B+H<sub>2</sub>]<sup>+</sup> ions ALWAYS follow their [M+H]<sup>+</sup> ions by a single MS scan, so a strong emphasis is placed on the RT scoring.   

The second component to our scoring system uses global modeling.  Global modeling serves as a method to ensure the key assumtions in the 1st component's scoring methods are correct. However, it also lets the researcher know about the overall quality of the group of scores produced for the putative DNA-adducts. The global model looks like this: 

</br></br>

<p align="center">
<img src="images/model.png" title = "scanning" width ="400">
</p>

</br></br></br>

<a name="getstart"/>


A simple likelihood value (like a probability) is produced from this statistical distribution of mass errors for the neutral loss.   The final score is simply the the probability multiplied by the previous feature score, very simple.


## Getting started

<a name="installation"/>

## 1. Installation
Start by downloading and installing the source R package. Don't forget to set your .libPaths() environment if needed:
```{r}

.libPaths("path to R library folder")

```
To finish installation, we strongly suggest downloading the repository as a zip package.   Then create a new R package using R-studio pointing to your un-zipped folder.   Then use the install button in RStudio to complete the installation into your R library.


## 2. Dependencies

#### Operating system:

Please note that due to the fact that software for acquiring and processing the mass spectrometry data used in this workflow was developed for Windows operating systems.   While some functions can be adapted to run on other operating systems, the very fact that Windows is used to run the mass spectrometers, to convert raw data to mzML,  and to process data into feature lists means you will likely need to run this workflow using a Windows computer.</br>

You can run a function that is provided in wSIMCity that will check for installed dependencies. These include:



#### R packages:

1. mzR : for reading mzML raw data
2. Rcpp : for mzR
3. doParallel : for multi threading
</br></br>


```{r}

wSIMCity::chckDependencies()

```
#### Windows compatible software:

1. The proteowizard msconvert software: for converting Thermo .Raw files (or any other supported vendor)
2. MSDIal : for feature finding.


## 3. Before you start<a name="beforestart"/>:

wSIMCity needs to have multiple items specified for it to work correctly.   

#### These include:

#### 1. A scan definition file.  

This file describes one duty cycle on the instrument in DIA SIM mode and defines what the *m/z* ranges for the wide SIM-MS<sup>2</sup> are.  It is tab delimited and is in the form:

|ScanType|WindowStart|WindowEnd|AquisitionStart|AcquisitionEnd|
|:---:|:---:|:---:|:---:|:---:|
|WSIM|197|364|330|364|
|NL|100|550|330|364|
|WSIM|197|394|360|394|
|NL|100|550|360|394|
|WSIM|197|424|390|424|
|NL|100|550|390|424|
|WSIM|197|454|420|454|
|NL|100|550|420|454|
|WSIM|197|484|450|484|
|NL|100|550|450|484|
|WSIM|197|514|480|514|
|NL|100|550|480|514|
|WSIM|197|544|510|544|
|NL|100|550|510|544|
|WSIM|197|574|540|574|
|NL|100|550|540|574|
|WSIM|197|604|570|604|
|NL|100|550|570|604|
|WSIM|197|634|600|634|
|NL|100|550|600|634|

##### Notes:

```ScanType``` is one of either 'WSIM' or 'NL' used to denote the scan level (MS<sup>1</sup> or MS<sup>2</sup>).</br>
```WindowStart``` and ```WindowEnd``` indicate the start and end *m/z* values for the data collection *m/z* range as set at the instrument.</br>
```AquisitionStart``` and ```AqcuisitionEnd``` denote the start and end *m/z* values for the mass range you filtered your data on during the run.

#### 2. A MSDIAL parameters file.

This file is the MSDIAL parameters file used with running the MSDIAL command line program.  wSIMCity will convert your raw data and then run MSDIAL to find features in your data.

#### 3. A list of adducts to search between MS scan levels.
You provide these as input to the neutral loss modelling step: modelNLM(adduct_mass = -116.0474), etc.


# Usage

## Set paths
```{r}

msconvert_path <- "C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.18229.34f38e1eb\\msconvert.exe"

msdial_path <- "C:\\MSTOOLS\\MSDIALv3.66\\MsdialConsoleApp.exe"

msdial_param_path <- "params\\msdial_param.txt"

results_dir <- "results"

raw_file_dir <- "raw"

scandef_file <- "params/scandef.txt"


```



## List your raw files and samples that are to be analyzed
```{r}

raw_file_list <- getRawFileList(raw_file_dir)

raw_file_list

sample_names <- getSampleNames(raw_file_dir)

sample_names


```


## prepare analysis results folders
```{r}

dir.create(results_dir)

makeSampleDir(results_dir,sample_names, scandef_file = scandef_file)

sample_directories <- getSampleDirectories(results_dir)

sample_directories



```




## convert raw file
```{r}

convertRaw(raw_file_dir = raw_file_dir, msconvert_path = msconvert_path)


```



## now segment the mzML datafile
```{r}

segmentMzMLDataSample(raw_file_dir, sample_directories, scandef_file = scandef_file)

```




## search mzml files with msdial
```{r}

findPeaksMSDIALSample(sample_directories,scandef_file,msdial_path,msdial_param_path,nCore = 1)

```


## Retrieve the MS-DIAL results
```{r}

msdial_results <- getMSDIAL_results(sample_directories)

save(file = "msdial_results.rda",msdial_results)


```


## Now perform the search on each file
A user friendly function to perform this loop will be available soon
```{r}
sample_directories


adduct_list <- data.frame("Neutral.Loss" = c("dR","13C_dR"),"MZ" = c(-116.0474,-121.0641))

system.time({
nlmMod <- modelNLM_run(msdial_results,
											 sample_directories,
											 adduct_list = adduct_list,
											 boost = 3,
											 alpha_mz = 0.5,
											 beta_rt = 0.5,
											 ppm_window = 30,
											 rt_tol = 0.3,
											 instrument_tol = .03)
})


```


## Filter and merge the results accross adducts and look at your table of results. 
```{r}


merge_results(sample_directories, adduct_list)


```



