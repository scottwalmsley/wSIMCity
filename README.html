<!doctype html><html><head><meta charset="utf-8">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/github-markdown-css/2.10.0/github-markdown.min.css">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/9.13.1/highlight.min.js">
<link  rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.10.0/dist/katex.min.css" integrity="sha384-9eLZqc9ds8eNjO3TmqPeYcDj8n+Qfa4nuSiGYa6DjLNcv9BtN69ZIulL9+8CqC9Y" crossorigin="anonymous">
<link rel="stylesheet" href="https://gitcdn.xyz/repo/goessner/mdmath/master/css/texmath.css">
<link rel="stylesheet" href="https://gitcdn.xyz/repo/goessner/mdmath/master/css/vscode-texmath.css">

</head><body class="markdown-body">
<h1 data-line="0" class="code-line" id="wsimcitybrbrsearching-for-dna-adducts-in-dia-wide-simmssup2sup-data">wSIMCity:</br></br>Searching for DNA-adducts in DIA wide-SIM/MS<sup>2</sup> data</h1>
</br>
<h2 data-line="4" class="code-line" id="table-of-contents">Table of Contents</h2>
<ul>
<li data-line="7" class="code-line">
<p data-line="7" class="code-line"><a href="#background">Background</a></p>
</li>
<li data-line="9" class="code-line">
<p data-line="9" class="code-line"><a href="#overview">Overview of software</a></br></p>
</li>
<li data-line="11" class="code-line">
<p data-line="11" class="code-line"><a href="#scoring">Scoring</a></p>
</li>
<li data-line="13" class="code-line">
<p data-line="13" class="code-line"><a href="#getstart">Getting started</a></p>
</li>
<li data-line="15" class="code-line">
<p data-line="15" class="code-line"><a href="#installation">Installation</a></br></p>
</li>
<li data-line="17" class="code-line">
<p data-line="17" class="code-line"><a href="#dependencies">Dependencies</a></br></p>
</li>
<li data-line="19" class="code-line">
<p data-line="19" class="code-line"><a href="#beforestart">Before you start</a></br></p>
</li>
</ul>
<p data-line="23" class="code-line"></br></br></p>
<h2 data-line="25" class="code-line" id="background">Background</h2>
<p data-line="30" class="code-line">Welcome to wSIMCity, software for mining Wide Sim mass spectrometry data.  wSIMCity was developed by Scott Walmsley, PhD, of the Masonic Cancer Center and the Institute for Health Informatics at the Univerity of Minnesota - Twin Cities for the laboratory of Prof. Rob Turesky, PhD.  Development of wSIMCity was brought about by necessity to overcome data structure issues in wide-SIM data independent acquisition (DIA) data produced in DNA-adduct mass spectrometry (MS) experiments, and to facilitate automated detection of DNA-adducts.  As such we developed an R-package to process raw wide-SIM DIA data and to mine for DNA-adducts in that data.</p>
<p data-line="32" class="code-line">Wide-SIM data independent acquisition (DIA) methods were developed to increase sensitivity for detection of DNA-adducts in digests of DNA and using mass spectrometry. More specfically, the automatic gain control (AGC) of orbitrap MS is used together with wide selected ion monitoring (wide-SIM) <em>m/z</em> ranges to boost sensitivity for the ions of interest. The wide-SIM MS<sup>1</sup> data is collected and then the subsequent scan undergoes a low energy collision induced dissociation (CID) at ~25 millielectron volts (meV).  This MS<sup>2</sup> scan leverages the known neutral loss at low energy CID of a dexoyribose and serves as the tell tale sign of the presence of a DNA-adduct.   Turesky and group pioneered this altered DIA methodology and it likely will serve as a strategy for the untargeted discovery of novel DNA-adducts. The general mechanism is shown below using the molecule dG-C8-PhIP, an important carcinogen of prostate cells caused by ingestion of overcooked meats:</br></p>
<p align="center">
<img src="images/mechanism.png" title = "mechanism">
</p>
<p data-line="39" class="code-line">The blue portion of the molecule is the mutagen or carcinogen, shown bound to the nucleoside (black and red).   MS<sup>2</sup> CID (25meV) causes the neutral loss of a deoxyribose (dR, red) leading to the formation of Guanine-C8-PhIP.  Generally speaking, it has been shown that this mechanism is universal across as wide spectrum of DNA adducts.   The ions at the MS<sup>1</sup> level are known as presursor ions (denoted: [M+H]<sup>+</sup>), and the resulting ion after the neutral loss is called the 'aglycone' (denoted as: [B+H<sub>2</sub>]<sup>+</sup>) because of the neutral loss of a deoxyribose.  Note the gain of a mobile proton on the aglycone molecule, hence the '+H<sub>2</sub>' notation.</p>
<p data-line="41" class="code-line">Targeted extraction of molecules such as dG-C8-PhIP from this type of data is fairly straightforward.   However, a major goal is to facilitate untargeted detection of DNA-adducts using this DIA based method. Below is a figure showing the slight differences between scanning methods common place in proteomics and metabolommics (A), and methods developed for DNA adductomics (B and C). Included in panel B is a DDA method for DNA adductomics  called constant neutral loss (CNL) screening.  Panel C is the method for which this software was originally written.  It is the very nature of the data structure produced in C that lead to incompatibilities with current DDA /DIA software developed for proteomics and metabolomics.</p>
<p data-line="43" class="code-line"></br></br></p>
<p align="center">
<img src="images/DIA_scanning.png" title = "scanning" width ="600">
</p>
<p data-line="49" class="code-line"></br></br></br></p>
<a name="overview"/>
<h2 data-line="53" class="code-line" id="overview-of-wsimcity">Overview of wSIMCity</h2>
<p data-line="55" class="code-line">wSIMCity seeks to discover the 'landscape' or map of DNA-adducts in a DNA sample prepared for and analyzed using wide-SIM-MS.  The overall workflow is shown in the next figure:</p>
<p data-line="57" class="code-line"></br></br></p>
<p align="center">
<img src="images/total_model.png" title = "scanning" width ="600">
</p>
<p data-line="64" class="code-line"></br></br></br></p>
<p data-line="67" class="code-line">Software developed for the DIA methods (and previously data dynamic acquisition (DDA) methodologies) have worked well, but are incompatible with our current wide-SIM scanning technique.   Therefor we have implemented a simple workaround using the R package <code>mzR</code>, incorporated the MSDIAL feature finding algortihm, and devised a custom global modelling and search strategy to detect the adducts.   Extracted ion chromatograms (EIC) for candidate molecules are produced.  Further we devised a scheme to best levereage a global scoring method together with an individual scoring metric to help the researcher identify and pursue candidate DNA-adduct molecules for identification.   The output of wSIMCity is list of masses and retention times to be used for targetted identification using high energy dissociation (HCD) MS.</p>
<h2 data-line="69" class="code-line" id="scoringa-name%22scoring%22">Scoring<a name="scoring"/></h2>
<p data-line="71" class="code-line">The best aspect of wSIMCIty is it produces a series of scores to help the researcher rummage through the list of candidate DNA-adducts.</p>
<p data-line="73" class="code-line">The score is broken down into several components. The first components are borrowed from MSDIAL's scoring method for matching detected compounds with those listed in a mass and retention database.  Instead, we use this scoreing method to match precursor - aglycone molecules found in DNA-adductomics data.</p>
<p data-line="75" class="code-line">A major difference between our scoring method and MSDIAL's is that the underlying assumptions about the distributive properties of the measured values in our data are a little different:  we use ppm mass errors to differentiate between real [M+H]<sup>+</sup> and [B+H<sub>2</sub>]<sup>+</sup> pairs of features and false hits. Specifically, for every [M+H]<sup>+</sup> ion, we compute a theroetical [M+H-dR]<sup>+</sup> <em>m/z</em>, and then look for [B+H<sub>2</sub>]<sup>+</sup> ions with that mass and retention time. We look for every ion with that <em>m/z</em> value and then measure how far off in ppm the experimental ion is. This measurement is called 'Mass err (ppm)' in the below plot.</p>
<p data-line="77" class="code-line">These measured data points follow a Laplace distribution (well, it's more like a Cauchy distribution, but that is near imposible to deal with in terms of the math involved). We also, like MSDIAL, incorporate retention time as a metric.  Like MSDIAL, our assumption here follows a Guassian distribution. [B+H<sub>2</sub>]<sup>+</sup> ions ALWAYS follow their [M+H]<sup>+</sup> ions by a single MS scan, so a strong emphasis is placed on the RT scoring.</p>
<p data-line="80" class="code-line">The first component scores can be summarized by the following equation:
</br></br></p>
<section><eqn><span class="katex-display"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msubsup><mi>S</mi><mrow><mi>N</mi><mi>L</mi></mrow><mrow><mi>u</mi><mi>n</mi><mi>k</mi></mrow></msubsup><mo>=</mo><mi>α</mi><mo>∗</mo><msup><mi>S</mi><msub><mi mathvariant="normal">Δ</mi><mrow><mi>p</mi><mi>p</mi><mi>m</mi></mrow></msub></msup><mo>+</mo><mi>β</mi><mo>∗</mo><msup><mi>S</mi><msub><mi mathvariant="normal">Δ</mi><mrow><mi>R</mi><mi>T</mi></mrow></msub></msup></mrow><annotation encoding="application/x-tex">S_{NL}^{unk} = \alpha * S^{\Delta_{ppm}} +  \beta * S^{\Delta_{RT}}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:1.146108em;vertical-align:-0.247em;"></span><span class="mord"><span class="mord mathdefault" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.8991079999999998em;"><span style="top:-2.4530000000000003em;margin-left:-0.05764em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.10903em;">N</span><span class="mord mathdefault mtight">L</span></span></span></span><span style="top:-3.113em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathdefault mtight">u</span><span class="mord mathdefault mtight">n</span><span class="mord mathdefault mtight" style="margin-right:0.03148em;">k</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.247em;"><span></span></span></span></span></span></span><span class="mspace" style="margin-right:0.2777777777777778em;"></span><span class="mrel">=</span><span class="mspace" style="margin-right:0.2777777777777778em;"></span></span><span class="base"><span class="strut" style="height:0.46528em;vertical-align:0em;"></span><span class="mord mathdefault" style="margin-right:0.0037em;">α</span><span class="mspace" style="margin-right:0.2222222222222222em;"></span><span class="mbin">∗</span><span class="mspace" style="margin-right:0.2222222222222222em;"></span></span><span class="base"><span class="strut" style="height:0.974661em;vertical-align:-0.08333em;"></span><span class="mord"><span class="mord mathdefault" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:0.891331em;"><span style="top:-3.113em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight"><span class="mord mtight">Δ</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.16454285714285716em;"><span style="top:-2.357em;margin-left:0em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.5em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight"><span class="mord mathdefault mtight">p</span><span class="mord mathdefault mtight">p</span><span class="mord mathdefault mtight">m</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.2818857142857143em;"><span></span></span></span></span></span></span></span></span></span></span></span></span></span></span><span class="mspace" style="margin-right:0.2222222222222222em;"></span><span class="mbin">+</span><span class="mspace" style="margin-right:0.2222222222222222em;"></span></span><span class="base"><span class="strut" style="height:0.8888799999999999em;vertical-align:-0.19444em;"></span><span class="mord mathdefault" style="margin-right:0.05278em;">β</span><span class="mspace" style="margin-right:0.2222222222222222em;"></span><span class="mbin">∗</span><span class="mspace" style="margin-right:0.2222222222222222em;"></span></span><span class="base"><span class="strut" style="height:0.8913309999999999em;vertical-align:0em;"></span><span class="mord"><span class="mord mathdefault" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:0.8913309999999999em;"><span style="top:-3.113em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight"><span class="mord mtight">Δ</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.3448em;"><span style="top:-2.3567071428571427em;margin-left:0em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.5em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.00773em;">R</span><span class="mord mathdefault mtight" style="margin-right:0.13889em;">T</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.14329285714285717em;"><span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></eqn></section><p data-line="85" class="code-line"></br></br></p>
<p data-line="87" class="code-line">....where ukn is the unknown feature, NL denotes 'neutral loss', ppm is the mass error in ppm, RT is retention time. <eq><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msub><mi mathvariant="normal">Δ</mi><mrow><mi>p</mi><mi>p</mi><mi>m</mi></mrow></msub></mrow><annotation encoding="application/x-tex">\Delta_{ppm}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.969438em;vertical-align:-0.286108em;"></span><span class="mord"><span class="mord">Δ</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.15139200000000003em;"><span style="top:-2.5500000000000003em;margin-left:0em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathdefault mtight">p</span><span class="mord mathdefault mtight">p</span><span class="mord mathdefault mtight">m</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.286108em;"><span></span></span></span></span></span></span></span></span></span></eq> is the measured error between the theoretical NL <em>m/z</em> value and the measured peak <em>m/z</em> in the MS<sup>2</sup> data.   <eq><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msub><mi mathvariant="normal">Δ</mi><mrow><mi>R</mi><mi>T</mi></mrow></msub></mrow><annotation encoding="application/x-tex">\Delta_{RT}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.83333em;vertical-align:-0.15em;"></span><span class="mord"><span class="mord">Δ</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.32833099999999993em;"><span style="top:-2.5500000000000003em;margin-left:0em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.00773em;">R</span><span class="mord mathdefault mtight" style="margin-right:0.13889em;">T</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.15em;"><span></span></span></span></span></span></span></span></span></span></eq>   is the measured difference between the MS<sup>1</sup> and MS<sup>2</sup> features. <eq><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>α</mi></mrow><annotation encoding="application/x-tex">\alpha</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.43056em;vertical-align:0em;"></span><span class="mord mathdefault" style="margin-right:0.0037em;">α</span></span></span></span></eq> and <eq><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>β</mi></mrow><annotation encoding="application/x-tex">\beta</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.8888799999999999em;vertical-align:-0.19444em;"></span><span class="mord mathdefault" style="margin-right:0.05278em;">β</span></span></span></span></eq> are values between 0-1 that weight the significance of the <eq><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msup><mi>S</mi><msub><mi mathvariant="normal">Δ</mi><mrow><mi>p</mi><mi>p</mi><mi>m</mi></mrow></msub></msup></mrow><annotation encoding="application/x-tex">S^{\Delta_{ppm}}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.841331em;vertical-align:0em;"></span><span class="mord"><span class="mord mathdefault" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:0.841331em;"><span style="top:-3.063em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight"><span class="mord mtight">Δ</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.16454285714285716em;"><span style="top:-2.357em;margin-left:0em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.5em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight"><span class="mord mathdefault mtight">p</span><span class="mord mathdefault mtight">p</span><span class="mord mathdefault mtight">m</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.2818857142857143em;"><span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></eq> or <eq><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msup><mi>S</mi><msub><mi mathvariant="normal">Δ</mi><mrow><mi>R</mi><mi>T</mi></mrow></msub></msup></mrow><annotation encoding="application/x-tex">S^{\Delta_{RT}}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.8413309999999999em;vertical-align:0em;"></span><span class="mord"><span class="mord mathdefault" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:0.8413309999999999em;"><span style="top:-3.063em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight"><span class="mord mtight">Δ</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.3448em;"><span style="top:-2.3567071428571427em;margin-left:0em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.5em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.00773em;">R</span><span class="mord mathdefault mtight" style="margin-right:0.13889em;">T</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.14329285714285717em;"><span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></eq> scores....they must add up to 1.  The independent scores are equated using:</p>
<p data-line="89" class="code-line"></br></br></p>
<section><eqn><span class="katex-display"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mstyle mathsize="1.44em"><msup><mi>S</mi><msub><mi mathvariant="normal">Δ</mi><mrow><mi>p</mi><mi>p</mi><mi>m</mi></mrow></msub></msup><mo>=</mo><mi>e</mi><mi>x</mi><msup><mi>p</mi><mrow><mo>−</mo><mo>(</mo><mfrac><mrow><mi mathvariant="normal">∣</mi><mi>o</mi><mi>b</mi><mi>s</mi><mo>−</mo><mi>r</mi><mi>e</mi><mi>f</mi><mi mathvariant="normal">∣</mi></mrow><mrow><mi>t</mi><mi>o</mi><mi>l</mi></mrow></mfrac><mo>)</mo></mrow></msup></mstyle></mrow><annotation encoding="application/x-tex">\Large S^{\Delta_{ppm}} = exp^{-(\frac{|obs-ref|}{tol})}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:1.27805em;vertical-align:0em;"></span><span class="mord sizing reset-size6 size8"><span class="mord mathdefault" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:0.8875347222222222em;"><span style="top:-3.413em;margin-right:0.034722222222222224em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size8 size6 mtight"><span class="mord mtight"><span class="mord mtight"><span class="mord mtight">Δ</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.15139200000000003em;"><span style="top:-2.5500000000000003em;margin-left:0em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathdefault mtight">p</span><span class="mord mathdefault mtight">p</span><span class="mord mathdefault mtight">m</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.286108em;"><span></span></span></span></span></span></span></span></span></span></span></span></span></span></span><span class="mspace" style="margin-right:0.2777777777777778em;"></span><span class="mrel sizing reset-size6 size8">=</span><span class="mspace" style="margin-right:0.2777777777777778em;"></span></span><span class="base"><span class="strut" style="height:1.8847135999999998em;vertical-align:-0.2799936em;"></span><span class="mord mathdefault sizing reset-size6 size8">e</span><span class="mord mathdefault sizing reset-size6 size8">x</span><span class="mord sizing reset-size6 size8"><span class="mord mathdefault">p</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:1.1143888888888887em;"><span style="top:-3.4130000000000003em;margin-right:0.034722222222222224em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size8 size6 mtight"><span class="mord mtight"><span class="mord mtight">−</span><span class="mopen mtight">(</span><span class="mord mtight"><span class="mopen nulldelimiter"></span><span class="mfrac"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:1.01em;"><span style="top:-2.6550000000000002em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathdefault mtight">t</span><span class="mord mathdefault mtight">o</span><span class="mord mathdefault mtight" style="margin-right:0.01968em;">l</span></span></span></span><span style="top:-3.23em;"><span class="pstrut" style="height:3em;"></span><span class="frac-line mtight" style="border-bottom-width:0.04em;"></span></span><span style="top:-3.485em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight">∣</span><span class="mord mathdefault mtight">o</span><span class="mord mathdefault mtight">b</span><span class="mord mathdefault mtight">s</span><span class="mbin mtight">−</span><span class="mord mathdefault mtight" style="margin-right:0.02778em;">r</span><span class="mord mathdefault mtight">e</span><span class="mord mathdefault mtight" style="margin-right:0.10764em;">f</span><span class="mord mtight">∣</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.345em;"><span></span></span></span></span></span><span class="mclose nulldelimiter"></span></span><span class="mclose mtight">)</span></span></span></span></span></span></span></span></span></span></span></span></span></eqn></section><p data-line="98" class="code-line"></br></br></p>
<section><eqn><span class="katex-display"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mstyle mathsize="1.44em"><msup><mi>S</mi><msub><mi mathvariant="normal">Δ</mi><mrow><mi>R</mi><mi>T</mi></mrow></msub></msup><mo>=</mo><mi>e</mi><mi>x</mi><msup><mi>p</mi><mrow><mo>−</mo><mfrac><mn>1</mn><mn>2</mn></mfrac><mo>(</mo><mfrac><mrow><mi mathvariant="normal">∣</mi><mi>M</mi><msubsup><mi>S</mi><mrow><mi>R</mi><mi>T</mi></mrow><mn>2</mn></msubsup><mo>−</mo><mi>M</mi><msubsup><mi>S</mi><mrow><mi>R</mi><mi>T</mi></mrow><mn>1</mn></msubsup><mi mathvariant="normal">∣</mi></mrow><mrow><mi>t</mi><mi>o</mi><mi>l</mi></mrow></mfrac><msup><mo>)</mo><mn>2</mn></msup></mrow></msup></mstyle></mrow><annotation encoding="application/x-tex">\Large S^{\Delta_{RT}} =  exp^{-\frac{1}{2}(\frac{|MS_{RT}^{2}-MS_{RT}^{1}|}{tol})^2}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:1.27805em;vertical-align:0em;"></span><span class="mord sizing reset-size6 size8"><span class="mord mathdefault" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:0.8875347222222222em;"><span style="top:-3.413em;margin-right:0.034722222222222224em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size8 size6 mtight"><span class="mord mtight"><span class="mord mtight"><span class="mord mtight">Δ</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.32833099999999993em;"><span style="top:-2.5500000000000003em;margin-left:0em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.00773em;">R</span><span class="mord mathdefault mtight" style="margin-right:0.13889em;">T</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.15em;"><span></span></span></span></span></span></span></span></span></span></span></span></span></span></span><span class="mspace" style="margin-right:0.2777777777777778em;"></span><span class="mrel sizing reset-size6 size8">=</span><span class="mspace" style="margin-right:0.2777777777777778em;"></span></span><span class="base"><span class="strut" style="height:2.2513526em;vertical-align:-0.2799936em;"></span><span class="mord mathdefault sizing reset-size6 size8">e</span><span class="mord mathdefault sizing reset-size6 size8">x</span><span class="mord sizing reset-size6 size8"><span class="mord mathdefault">p</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:1.3689993055555556em;"><span style="top:-3.4130000000000003em;margin-right:0.034722222222222224em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size8 size6 mtight"><span class="mord mtight"><span class="mord mtight">−</span><span class="mord mtight"><span class="mopen nulldelimiter"></span><span class="mfrac"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.845108em;"><span style="top:-2.6550000000000002em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight">2</span></span></span></span><span style="top:-3.23em;"><span class="pstrut" style="height:3em;"></span><span class="frac-line mtight" style="border-bottom-width:0.04em;"></span></span><span style="top:-3.394em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight">1</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.345em;"><span></span></span></span></span></span><span class="mclose nulldelimiter"></span></span><span class="mopen mtight">(</span><span class="mord mtight"><span class="mopen nulldelimiter"></span><span class="mfrac"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:1.376639em;"><span style="top:-2.6550000000000002em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathdefault mtight">t</span><span class="mord mathdefault mtight">o</span><span class="mord mathdefault mtight" style="margin-right:0.01968em;">l</span></span></span></span><span style="top:-3.23em;"><span class="pstrut" style="height:3em;"></span><span class="frac-line mtight" style="border-bottom-width:0.04em;"></span></span><span style="top:-3.623831em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight">∣</span><span class="mord mathdefault mtight" style="margin-right:0.10903em;">M</span><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:1.07544em;"><span style="top:-2.2516700000000003em;margin-left:-0.05764em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.7em;"></span><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.00773em;">R</span><span class="mord mathdefault mtight" style="margin-right:0.13889em;">T</span></span></span><span style="top:-3.1310000000000002em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.7em;"></span><span class="mord mtight"><span class="mord mtight">2</span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.44833em;"><span></span></span></span></span></span></span><span class="mbin mtight">−</span><span class="mord mathdefault mtight" style="margin-right:0.10903em;">M</span><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.05764em;">S</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:1.07544em;"><span style="top:-2.2516700000000003em;margin-left:-0.05764em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.7em;"></span><span class="mord mtight"><span class="mord mathdefault mtight" style="margin-right:0.00773em;">R</span><span class="mord mathdefault mtight" style="margin-right:0.13889em;">T</span></span></span><span style="top:-3.1310000000000002em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.7em;"></span><span class="mord mtight"><span class="mord mtight">1</span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.44833em;"><span></span></span></span></span></span></span><span class="mord mtight">∣</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.345em;"><span></span></span></span></span></span><span class="mclose nulldelimiter"></span></span><span class="mclose mtight"><span class="mclose mtight">)</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:0.8141079999999999em;"><span style="top:-3.063em;margin-right:0.05em;"><span class="pstrut" style="height:2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight">2</span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></eqn></section><p data-line="102" class="code-line"></br></br>
...where <eq><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>t</mi><mi>o</mi><mi>l</mi></mrow><annotation encoding="application/x-tex">tol</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:0.69444em;vertical-align:0em;"></span><span class="mord mathdefault">t</span><span class="mord mathdefault">o</span><span class="mord mathdefault" style="margin-right:0.01968em;">l</span></span></span></span></eq> is a user defined instrument mass tolerance window in ppm (<em>e.g.</em> 5 ppm) or expected max deviation of retention times for the MS<sup>1</sup> and MS<sup>2</sup> features.</p>
<p data-line="107" class="code-line">The second component to our scoring system uses global modeling.  Global modeling serves as a method to ensure the key assumtions in the 1st component's scoring methods are correct. However, it also lets the researcher know about the overall quality of the group of scores produced for the putative DNA-adducts. The global model looks like this:</p>
<p data-line="109" class="code-line"></br></br></p>
<p align="center">
<img src="images/model.png" title = "scanning" width ="400">
</p>
<p data-line="115" class="code-line"></br></br></br></p>
<a name="getstart"/>
<h2 data-line="119" class="code-line" id="getting-started">Getting started</h2>
<a name="installation"/>
<h2 data-line="123" class="code-line" id="1-installation">1. Installation</h2>
<p data-line="124" class="code-line">Start by downloading and installing the source R package. Don't forget to set your .libPaths() environment if needed.</p>
<pre><code data-line="125" class="code-line language-{r}"><code><div>
.libPaths(&quot;path to R library folder&quot;)

devtools::install_github(&quot;scottwalmsley/wSIMCity&quot;)


</div></code></code></pre>
<h2 data-line="134" class="code-line" id="2-dependencies">2. Dependencies</h2>
<h4 data-line="136" class="code-line" id="operating-system">Operating system:</h4>
<p data-line="138" class="code-line">Please note that due to the fact that software for acquiring and processing the mass spectrometry data used in this workflow was developed for Windows operating systems.   While some functions can be adapted to run on other operating systems, the very fact that Windows is used to run the mass spectrometers, to convert raw data to mzML,  and to process data into feature lists means you will likely need to run this workflow using a Windows computer.</br></p>
<p data-line="140" class="code-line">You can run a function that is provided in wSIMCity that will check for installed dependencies. These include:</p>
<h4 data-line="144" class="code-line" id="r-packages">R packages:</h4>
<ol>
<li data-line="146" class="code-line">mzR : for reading mzML raw data</li>
<li data-line="147" class="code-line">Rcpp : for mzR</li>
<li data-line="148" class="code-line">doParallel : for multi threading
</br></br></li>
</ol>
<pre><code data-line="152" class="code-line language-{r}"><code><div>
wSIMCity::chckDependencies()

</div></code></code></pre>
<h4 data-line="157" class="code-line" id="windows-compatible-software">Windows compatible software:</h4>
<ol>
<li data-line="159" class="code-line">The proteowizard msconvert software: for converting Thermo .Raw files (or any other supported vendor)</li>
<li data-line="160" class="code-line">MSDIal : for feature finding.</li>
</ol>
<h2 data-line="163" class="code-line" id="3-before-you-starta-name%22beforestart%22">3. Before you start<a name="beforestart"/>:</h2>
<p data-line="165" class="code-line">wSIMCity needs to have multiple items specified for it to work correctly.</p>
<h4 data-line="167" class="code-line" id="these-include">These include:</h4>
<h4 data-line="169" class="code-line" id="1-a-scan-definition-file">1. A scan definition file.</h4>
<p data-line="171" class="code-line">This file describes one duty cycle on the instrument in DIA SIM mode and defines what the <em>m/z</em> ranges for the wide SIM-MS<sup>2</sup> are.  It is tab delimited and is in the form:</p>
<table>
<thead>
<tr>
<th style="text-align:center">ScanType</th>
<th style="text-align:center">WindowStart</th>
<th style="text-align:center">WindowEnd</th>
<th style="text-align:center">AquisitionStart</th>
<th style="text-align:center">AcquisitionEnd</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">364</td>
<td style="text-align:center">330</td>
<td style="text-align:center">364</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">330</td>
<td style="text-align:center">364</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">394</td>
<td style="text-align:center">360</td>
<td style="text-align:center">394</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">360</td>
<td style="text-align:center">394</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">424</td>
<td style="text-align:center">390</td>
<td style="text-align:center">424</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">390</td>
<td style="text-align:center">424</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">454</td>
<td style="text-align:center">420</td>
<td style="text-align:center">454</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">420</td>
<td style="text-align:center">454</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">484</td>
<td style="text-align:center">450</td>
<td style="text-align:center">484</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">450</td>
<td style="text-align:center">484</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">514</td>
<td style="text-align:center">480</td>
<td style="text-align:center">514</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">480</td>
<td style="text-align:center">514</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">544</td>
<td style="text-align:center">510</td>
<td style="text-align:center">544</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">510</td>
<td style="text-align:center">544</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">574</td>
<td style="text-align:center">540</td>
<td style="text-align:center">574</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">540</td>
<td style="text-align:center">574</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">604</td>
<td style="text-align:center">570</td>
<td style="text-align:center">604</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">570</td>
<td style="text-align:center">604</td>
</tr>
<tr>
<td style="text-align:center">WSIM</td>
<td style="text-align:center">197</td>
<td style="text-align:center">634</td>
<td style="text-align:center">600</td>
<td style="text-align:center">634</td>
</tr>
<tr>
<td style="text-align:center">NL</td>
<td style="text-align:center">100</td>
<td style="text-align:center">550</td>
<td style="text-align:center">600</td>
<td style="text-align:center">634</td>
</tr>
</tbody>
</table>
<h5 data-line="196" class="code-line" id="notes">Notes:</h5>
<p data-line="198" class="code-line"><code>ScanType</code> is one of either 'WSIM' or 'NL' used to denote the scan level (MS<sup>1</sup> or MS<sup>2</sup>).</br>
<code>WindowStart</code> and <code>WindowEnd</code> indicate the start and end <em>m/z</em> values for the data collection <em>m/z</em> range as set at the instrument.</br>
<code>AquisitionStart</code> and <code>AqcuisitionEnd</code> denote the start and end <em>m/z</em> values for the mass range you filtered your data on during the run.</p>
<h4 data-line="202" class="code-line" id="2-a-msdial-parameters-file">2. A MSDIAL parameters file.</h4>
<p data-line="204" class="code-line">This file is the MSDIAL parameters file used with running the MSDIAL command line program.  wSIMCity will convert your raw data and then run MSDIAL to find features in your data.</p>
<h4 data-line="206" class="code-line" id="3-a-file-of-adducts-to-search-between-ms-scan-levels">3. A file of adducts to search between MS scan levels.</h4>
<p data-line="208" class="code-line">This is the tab delimited text file containing the list of adducts masses you wish to search.</p>
<table>
<thead>
<tr>
<th style="text-align:center">Neutral Loss</th>
<th style="text-align:center">MZ</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center">dR</td>
<td style="text-align:center">-116.0474</td>
</tr>
<tr>
<td style="text-align:center">[<sup>13</sup>C]-dR</td>
<td style="text-align:center">-121.0641</td>
</tr>
</tbody>
</table>
<h5 data-line="216" class="code-line" id="set-your-paths-and-environment-variables">Set your paths and environment variables.</h5>
<pre><code data-line="218" class="code-line language-{r}"><code><div>



</div></code></code></pre>

</body></html>