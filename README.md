# BulkSeqAssistant V5.2  
by Zhiming Ye @ Yao Lab, SKLRD   
  
> Part of the code comes from online materials. For learning and communication purposes only.  


### A function collection developed to simplify the omics analysis process.
This package is designed for bulk RNA-seq analysis and organizing data.  
Part of the function can be adopt in organizing data, WES analysis and scRNA-seq analysis as well.  
### NEWs in V5.2 (2022-09-12):  
- Add a function : BuildMultigroupDEGlist.  
- The orgDB parameter enable user to select the annotation database.  
- Fix bugs in PlotNES.  
- 'CorrlatinEstimate0' function is now replaced by 'CorrEstimate0'
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
