---
title: "README"
author: "Scott Walmsley"
date: "May 3, 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## wSIMCity: mining for DNA-adducts using wide - SIM mass spectrometry

Welcome to wSIMCity, software for mining Wide Sim mass spectrometry data.  wSIMCity was developed by Scott Walmsley, PhD, of the Masonic Cancer Center at the Univerity of Minnesota - Twin Cities in the laboratory of Prof. Rob Turesky, PhD.  Development of wSIMCity was brought about by necessity to overcome data structure issues in wide-SIM data independent acquisition (DIA) data produced in DNA-adduct mass spectrometry (MS) experiments.  As such we developed an R-package to process raw w-SIM DIA data and to mine for DNA-adducts in that data.   wSIM-DIA methods were developed to increase sensitivities 

## 1. Installation
Start by downloading and installing the source R package.
```{r}

devtools::install_github("scottwalmsley/wSIMCity")


```




