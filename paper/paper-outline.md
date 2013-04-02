#  Paper Outline

## INTRODUCTION

### Microarray Design & Work flow 

* Intro to basic idea of 450k testing for DNA methylation
 
 * Isolation, preparation, bisulphite conversion
 * 2 probes for each CpG
 * What is a beta value
  
* Introduce type I probes and 27k chip
 * Single base extension on 2 probe types
 * 1 color channel
 
* Introduce 450k & type II probes and compare/contrast
 * Competitive hybridization & extension for single probe type
 * 2 color channels
 
* Red and green channels for type I and type II **table for 6 signal types**
 * Which color channels are added for which base
 * Signal characteristics/qualities **horseshoe plots**

### Motivation for milr normalization

 * Cross hybridization, etc
 * Multiple probe types with very different characteristics
 * Cell type mixture (i.e. brain tissue)

### Introduction of normalization method 	
* Local regression for intensity dependence on BOTH channels & GC content
 * We model deviations from the experiment wise average as function of these
 * Insert Y = Ax +b type **formula** here for channels  
 
* House keeping controls are invariant to cell mixtures
 * Critical assumption: most HK CpGs are similar between samples
  * May not be true for some cancers
  * Talk about HK CpG selection, maybe include figure or leave for a supplement

## METHODS

*	Math & Stats part
 * Steps of algorithm: take log2, median center, etc…
 * Brief description of span selection & local regression surface **formula**
 
*	Tedious stratification part
 * How do we split apart and normalize the different signals **flow chart**

## RESULTS
*	Introduce data sets
 * BrainSpan
 * CFS 
 
* Present results
 * Before & after clustering **plots using milr, SWAN, and IMA peak correction**
 * Before & after p-values for batch effect & dependent variable diffs **plot/table**

## DISCUSSION






