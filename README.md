# BulkSeqAssistant V5.3  
by Zhiming Ye @ Yao Lab, SKLRD   
  
> Part of the code comes from online materials. For learning and communication purposes only.  


### A function collection developed to simplify the omics analysis process.
This package is designed for bulk RNA-seq analysis and organizing data.  
Part of the function can be adopt in organizing data, WES analysis and scRNA-seq analysis as well.  
### NEWs in V5.3 (2022-10-02):  
- Add a function : MultiAOV, for Multi-group parametric test.  
### Installation
```
library(devtools)
install_github("ZhimingYe/BulkSeqAssistant")
```
WARNING: BioConductor dependencies like 'clusterProfiler', 'edgeR', 'GSEABase', 'GSVA', 'ReactomePA', 'enrichplot', 'impute' and 'preprocessCore', as well as 'ggbiplot' form Github were not on CRAN. So you may install them by yourself.  
### Function List: 
* Arrange_Table
* BuildMultigroupDEGlist
* CalcMad
* CorrEstimate0
* do.GSVA
* do.ssGSEA
* doGO
* doKEGG
* doRA
* fpkmToTpm
* getDESeq2DEG
* MapSymbol
* MultiDunn
* NMFEstimate
* NormalizeTest
* Plot.Box.Lite
* Plot.lm.smooth
* PlotNES
* PlotNES.lite
* PlotPCA
* PlotPie
* PlotSanky
* PlotUMAP
* PrettyMsigdbName
* Print.Char
* Ranked.GS
* rmEnsemblDot
* rmSep
* UniCox
* WGCNACalcSFT
* wilcox.adj.test
* wilDEGs
* Write.tables
