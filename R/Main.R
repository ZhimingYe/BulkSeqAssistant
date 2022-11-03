#Last changed : 2022-10-02 0:14

#library(tidyverse)
#library(magrittr)
#' @title Get DEG from DEseq2 result
#'
#' @param DEseq2Result A DEseq2 result. Should be converted to data frame.
#' @param Type Can be "all", "up" and "down".
#' @param Padj P adjusted cut off value.
#' @param lFC logFoldChange cut off value.
#' @param limma DEseq2 result is False, but limma result please enter T.
#'
#' @return A filtered DEG data frame.
#' @export
#' @author Zhiming Ye
#'
#' @examples
getDESeq2DEG<-function(DEseq2Result,Type,Padj=0.05,lFC=1,limma=F){
  if(limma){
    DEseq2Result<-DEseq2Result%>%dplyr::rename(log2FoldChange=logFC,padj=adj.P.Val)
  }
  DEseq2Result0<-DEseq2Result%>%dplyr::filter(padj<Padj)%>%dplyr::filter(abs(log2FoldChange)>lFC)
  if(Type=="all"){
    return(DEseq2Result0)
  }else if(Type=="up"){
    DEseq2Result0<-DEseq2Result0%>%dplyr::filter(log2FoldChange>lFC)
  }else if(Type=="down"){
    DEseq2Result0<-DEseq2Result0%>%dplyr::filter(log2FoldChange<(-lFC))
  }else{}
}

#' @title An enhanced method to print character vector to the console
#'
#' @param CharacterCollection A character vector which will be print.
#' @param Type A character. Can be "c", "tab" or "plus". c means comma, tab means print in multi-lines, and plus means to separate with "+"
#'
#' @return
#' @export
#' @author Zhiming Ye
#'
#' @examples
Print.Char<-function(CharacterCollection,Type="c"){

  if(Type=="tab"){
    CHAR<-"\n"
    for(i in 1:length(CharacterCollection)){
      if(i!=length(CharacterCollection)){
        CHAR<-paste0(CHAR,CharacterCollection[i],"\n")
      }
      else{
        CHAR<-paste0(CHAR,CharacterCollection[i])
      }
    }
  }
  if(Type=="c"){
    CHAR<-"\""
    for(i in 1:length(CharacterCollection)){
      if(i!=length(CharacterCollection)){
        CHAR<-paste0(CHAR,CharacterCollection[i],"\",\"")
      }
      else{
        CHAR<-paste0(CHAR,CharacterCollection[i],"\"")
      }
    }
  }
  if(Type=="plus"){
    CHAR<-"~"
    for(i in 1:length(CharacterCollection)){
      if(i!=length(CharacterCollection)){
        CHAR<-paste0(CHAR,CharacterCollection[i],"+")
      }
      else{
        CHAR<-paste0(CHAR,CharacterCollection[i])
      }
    }
  }
  cat(CHAR)
}

#' @title Two-group Box plot with stastical test
#'
#' @param Mat A table
#' @param xcolName Name of group 1
#' @param ycolName Name of group 2
#' @param colorList A vector of colors
#' @param my_comparisons Comparison list
#' @param Method like "wilcox.test" and "t.test"
#'
#' @return A ggplot2 object
#' @export
#' @author Zhiming Ye
#'
#' @examples
Plot.Box.Lite<-function(Mat,xcolName,ycolName,colorList,my_comparisons,Method="wilcox.test"){
  p<-ggboxplot(
    Mat,
    x = xcolName,
    y = ycolName,
    fill = xcolName,
    bxp.errorbar = T,
    palette = colorList
  ) + guides(fill = "none") + stat_compare_means(comparisons = my_comparisons,
                                                 method = Method,
                                                 label = "p.signif") + rotate_x_text(45)
  return(p)
}

#' @title A DEG analysis method by wilcox test
#' @description input1:data_matrix，rowname is gene，colname is sample; input2:sample_group(df),column must contains Sample_ID and Type
#' @param expr_matrix
#' @param sample_df
#' @param pval
#' @param fdr
#' @param FC
#' @param adj.method
#'
#' @return
#' @export
#' @author Zhiming Ye
#'
#' @examples
wilDEGs <- function(expr_matrix, sample_df, pval, fdr, FC, adj.method = "BH"){
  # tot_sam_ct <- ncol(expr_matrix)
  # tot_gene_ct <- nrow(expr_matrix)
  # tot_sam_label <- nrow(sample_df)
  # cat("input1:data_matrix，rowname is gene，colname is sample\ninput2:sample_group(df),column must contains Sample_ID and Type\n")
  # 提取shared samples
  if("Sample_ID" %in% colnames(sample_df)){
    inter_sam <- intersect(colnames(expr_matrix),
                           sample_df$Sample_ID)
  }else{
    return("Sample_df should contain two cols: Sample_ID,Type")
  }
  if(length(inter_sam) == 0){
    return("No samples shared in expr data and sample_df")
  }else{
    ana_expr <- expr_matrix[,as.character(inter_sam)]
    ana_df <- sample_df %>% filter(Sample_ID %in% inter_sam)
    cat("Shared sample count ",nrow(ana_df),"\n")
  }

  # 定义输出结果
  out_data <- data.frame(Gene = rownames(ana_expr),
                         MeanA = 0,MeanB = 0,
                         LogFC = 0,Pvalue = 1,FDR = 1)
  type1 <- unique(ana_df$Type)[1]
  type2 <- unique(ana_df$Type)[2]

  colnames(out_data)[2] <- paste0("Mean_",type1)
  colnames(out_data)[3] <- paste0("Mean_",type2)

  # 逐列进行计算均值，FC，P值，FDR
  out_data[,2] <- apply(ana_expr,1,function(x){
    mean(x[ana_df$Sample_ID[ana_df$Type == type1]],na.rm = T)
  })
  cat("Finished cal of mean ",type1,"\n")
  out_data[,3] <- apply(ana_expr,1,function(x){
    mean(x[ana_df$Sample_ID[ana_df$Type == type2]],na.rm = T)
  })
  cat("Finished cal of mean ",type2,"\n")
  out_data$Pvalue <- apply(ana_expr,1,function(x){
    wl <- wilcox.test(x[ana_df$Sample_ID[ana_df$Type == type1]],
                      x[ana_df$Sample_ID[ana_df$Type == type2]])
    # round(wl$p.value,4)
    wl$p.value
  })
  cat("Finished cal of Pvalue\n")
  out_data$FDR <- p.adjust(out_data$Pvalue, method = adj.method)

  out_data$LogFC <- round(log2(out_data[,3]/out_data[,2]),3)

  # 按照阈值生成最后的label
  out_data <- out_data %>%
    mutate(Label = ifelse(
      LogFC >= log2(FC) & Pvalue <= pval &
        FDR <= fdr,paste0(type2,"_up"),"Stable")) %>%
    mutate(Label = ifelse(
      LogFC <= -log2(FC) & Pvalue <= pval &
        FDR <= fdr,paste0(type2,"_down"),Label))

  return(out_data)
}

#' @title Plot the fitted line
#'
#' @param dataset variables should be columns.
#' @param xval
#' @param yval
#' @param Method like "pearson" or "spearman"
#'
#' @return A ggplot2 object
#' @export
#' @author Zhiming Ye
#'
#' @examples
Plot.lm.smooth<-function(dataset,xval,yval,Method="spearman"){
  #library(ggplot2)
  #library(ggpubr)
  p<-ggplot(data=dataset, aes(x=!!sym(xval), y=!!sym(yval)))+geom_point(color="black")+stat_smooth(method="lm",se=T)+stat_cor(data=dataset, method = Method)+theme_bw()
  return(p)
}
#' @title Return ranked gene list which is use for "GSEA" and "cnetplot" in package "clusterProfiler"
#'
#' @param Gene A vector containing genes
#' @param log2FC A vector containg log2FC
#' @param FromType default is "SYMBOL"
#' @param OrgDB default as org.Hs.eg.db
#'
#' @return A ranked named numeric vector. Names of the numbers is the ENTREZID.
#' @export
#' @author Zhiming Ye
#'
#' @examples
Ranked.GS<-function(Gene,log2FC,FromType = "SYMBOL",OrgDB=org.Hs.eg.db){
  Genetable<-data.frame(Gene=Gene,log2FC=log2FC)
  #library(clusterProfiler)
  #library(enrichplot)
  #library(ReactomePA)
  #library(org.Hs.eg.db)
  ENTREZIDtable<-clusterProfiler::bitr(Genetable$Gene,fromType = FromType,toType = "ENTREZID",OrgDb = OrgDB)
  colnames(ENTREZIDtable)[1]<-"Gene"
  Genetable<-Genetable%>%left_join(ENTREZIDtable)%>%arrange(desc(log2FC))
  GSElist<-as.numeric(Genetable$log2FC)
  names(GSElist)<-Genetable$ENTREZID
  GSElist = sort(GSElist, decreasing = TRUE)
  return(GSElist)
}
#' @title ORA GO analysis by clusterProfiler by symbol identifiers.
#'
#' @param GS Gene Symbols for enrichment analysis. NOT ENSEMBL ID.
#' @param ont One of "BP", "CC", "MF"  and "ALL".
#' @param PVal P adjusted value cut off, default 0.01
#' @param QVal qvalue cut off, default 0.05
#' @param OrgDB default as org.Hs.eg.db
#'
#' @return An enrichment result generated by clusterProfiler.
#' @export
#' @author Zhiming Ye
#'
#' @examples
doGO<-function(GS,ont="BP",PVal=0.01,QVal=0.05,OrgDB=org.Hs.eg.db){
  #library(clusterProfiler)
  #library(enrichplot)
  #library(ReactomePA)
  #library(org.Hs.eg.db)
  GeneList<-GS%>%unique()%>%bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = OrgDB)
  EnrichGO<-enrichGO(GeneList$ENTREZID,ont = ont,OrgDb = OrgDB,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable(OrgDb=OrgDB,keyType = "ENTREZID")
  return(EnrichGO)
}
#' @title ORA KEGG analysis by clusterProfiler by symbol identifiers.
#'
#' @param GS Gene Symbols for enrichment analysis. NOT ENSEMBL ID.
#' @param PVal P adjusted value cut off, default 0.01
#' @param QVal qvalue cut off, default 0.05
#' @param OrgDB default as org.Hs.eg.db
#'
#' @return An enrichment result generated by clusterProfiler.
#' @export
#' @author Zhiming Ye
#'
#' @examples
doKEGG<-function(GS,PVal=0.01,QVal=0.05,OrgDB=org.Hs.eg.db){
  #library(clusterProfiler)
  #library(enrichplot)
  #library(ReactomePA)
  #library(org.Hs.eg.db)
  GeneList<-GS%>%unique()%>%bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = OrgDB)
  EnrichKEGG<-enrichKEGG(GeneList$ENTREZID,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable(OrgDb=OrgDB,keyType = "ENTREZID")
  return(EnrichKEGG)
}
#' @title ORA Reactome analysis by clusterProfiler by symbol identifiers.
#'
#' @param GS Gene Symbols for enrichment analysis. NOT ENSEMBL ID.
#' @param PVal P adjusted value cut off, default 0.01
#' @param QVal qvalue cut off, default 0.05
#' @param OrgDB default as org.Hs.eg.db
#'
#' @return An enrichment result generated by clusterProfiler.
#' @export
#' @author Zhiming Ye
#'
#' @examples
doRA<-function(GS,PVal=0.01,QVal=0.05,OrgDB=org.Hs.eg.db){
  #library(clusterProfiler)
  #library(enrichplot)
  #library(ReactomePA)
  #library(org.Hs.eg.db)
  GeneList<-GS%>%unique()%>%bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = OrgDB)
  EnrichRA<-enrichPathway(GeneList$ENTREZID,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable(OrgDb=OrgDB,keyType = "ENTREZID")
  return(EnrichRA)
}

#' @title A simple PCA plot generator.
#'
#' @param Mat Expression levels.
#' @param group A vector, The order shall be consistent with the table.
#' @param Ellipse
#' @param counts If the expression table is original counts, should be TRUE. If is the TPM, FPKM or RPKM, should be FALSE.
#' @param loged Whether the expression level is log-transformed.
#' @param NumOfGenes Choose the top n genes with the highest MAD value for PCA analysis.
#'
#' @return A ggplot2 object
#' @export
#' @author Zhiming Ye
#'
#' @examples
PlotPCA<-function(Mat,group,Ellipse=T,counts=F,loged=T,NumOfGenes=10000){
  #library(FactoMineR)
  #library(factoextra)
  #library(ggbiplot)
  #library(edgeR)
  if(counts){
    Mat<-log2(edgeR::cpm(Mat)+1)
  }
  if((!loged)&(!counts)){
    Mat<-log2(Mat+1)
  }
  Mat<-CalcMad(Mat,NumOfGenes)
  res.pca2 <- prcomp(Mat%>%t(),scale=T)
  pic1<-ggbiplot.internal(res.pca2,var.axes = F,obs.scale = 0.5,groups = as.factor(group),ellipse = Ellipse,circle = F)+theme_bw()
  return(pic1)
}

#' @title A simple UMAP plot generator.
#'
#' @param Mat Expression levels.
#' @param group A vector, The order shall be consistent with the table.
#' @param counts If the expression table is original counts, should be TRUE. If is the TPM, FPKM or RPKM, should be FALSE.
#' @param loged Whether the expression level is log-transformed.
#' @param PointSize Point size of figure.
#' @param NumOfGenes Choose the top n genes with the highest MAD value for PCA analysis.
#'
#' @return A ggplot2 object
#' @export
#' @author Zhiming Ye
#'
#' @examples
PlotUMAP<-function(Mat,group,counts=F,loged=T,PointSize=2,NumOfGenes=10000){
  #library(umap)
  #library(edgeR)
  if(counts){
    Mat<-log2(edgeR::cpm(Mat)+1)
  }
  if((!loged)&(!counts)){
    Mat<-log2(Mat+1)
  }
  Mat<-CalcMad(Mat,NumOfGenes)
  res.umap3 <- umap::umap(Mat%>%t())
  draw.UMAP<-res.umap3[["layout"]]%>%as.data.frame()%>%dplyr::mutate(Group=group)
  colnames(draw.UMAP)[c(1,2)]<-c("UMAP_1","UMAP_2")
  pic<-ggplot(draw.UMAP,aes(UMAP_1,UMAP_2,color=Group))+geom_point(size=PointSize)+theme_bw()+labs(title = "UMAP visualization")
  return(pic)
}

#' @title Batch wilCox test, with P value Corrected.
#'
#' @param Mat Columns should be every variable to be tested, and another one column contains information about grouping, which is refer to the FactorCol
#' @param FactorCol Which column contains the grouping information. The number of it.
#' @param Method Default is "BH".
#'
#' @return
#' @export
#' @author Zhiming Ye
#'
#' @examples
wilcox.adj.test<-function(Mat,FactorCol,Method="BH"){
  WillTest<-Mat
  ResultRes<-data.frame(PValue=c(1),CellType=c("CT"))
  colnames(WillTest)[FactorCol]<-"Cluster"
  CellTypeList<-colnames(WillTest)[-FactorCol]
  for(CellType in CellTypeList){
    TestSet<-WillTest[,c("Cluster",CellType)]
    colnames(TestSet)[2]<-"Cell"
    KwRes<-wilcox.test(Cell~Cluster,TestSet)
    ResDf0<-data.frame(PValue=KwRes[["p.value"]],CellType=CellType)
    ResultRes<-rbind(ResultRes,ResDf0)
  }
  ResultRes<-ResultRes[-1,]
  ResultRes<-ResultRes%>%dplyr::mutate(Padj=p.adjust(.$PValue,method=Method))
  return(ResultRes)
}

#' @title Multi group nonparametric test
#'
#' @param Mat Columns should be every variable to be tested, and another one column contains information about grouping, which is refer to the FactorCol
#' @param FactorCol Which column contains the grouping information. The number of it.
#' @param Method Default the "bonferroni"
#' @param OnlySig Whether show only significant result.
#' @param pcutoff P value cut off.
#'
#' @return
#' @export
#' @author Zhiming Ye
#'
#' @examples
MultiDunn<-function(Mat,FactorCol,Method="bonferroni",OnlySig=F,pcutoff=0.01){
  #library(FSA)
  WillTest<-Mat
  ResultRes<-data.frame()
  colnames(WillTest)[FactorCol]<-"Cluster"
  CellTypeList<-colnames(WillTest)[-FactorCol]
  for(CellType in CellTypeList){
    TestSet<-WillTest[,c("Cluster",CellType)]
    colnames(TestSet)[2]<-"Cell"
    KwRes<-kruskal.test(Cell~Cluster,TestSet)

    if(KwRes[["p.value"]]<pcutoff){
      cat(paste0("kruskal test PASS : ",CellType," , Processing Dunn Test...\n"))
      WRes<-dunnTest(Cell~Cluster,TestSet,method=Method)
      Result0<-WRes[["res"]]%>%dplyr::mutate(GroupName=CellType)
      ResultRes<-rbind(ResultRes,Result0)
    }
  }
  if(OnlySig){
    ResultRes2<-ResultRes%>%dplyr::filter(P.adj<pcutoff)
    ResultRes3<-ResultRes%>%dplyr::filter(GroupName%in%names(table(ResultRes2$GroupName)))
    return(ResultRes3)
  }
  else{
    return(ResultRes)
  }
}

#' @title Generate a GSEA plot with ranked curve, which only used when drawing multi-pathways.
#'
#' @param gseRes GSEA result generated by clusterProfiler
#' @param GeneSet Description of the selected pathway. A vector.
#' @param nCols layout of figure, numbers of columns.
#' @param AnnoTerms Default the NES and p.adjust
#'
#' @return A cowplot plot_gird object
#' @export
#' @author Zhiming Ye
#'
#' @examples
PlotNES.lite<-function(gseRes,GeneSet,nCols=2,AnnoTerms= c("NES","p.adjust")){
  pp <- lapply(which(gseRes@result$Description%in%GeneSet), function(i) {
    anno <- gseRes[i,AnnoTerms]
    lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

    gsearank(gseRes, i, gseRes[i, 2]) + xlab(NULL) +ylab(NULL) +
      annotate("text", 10000, gseRes[i, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)
  })
  pp2<-plot_grid(plotlist=pp, ncol=nCols,align = "hv")
  return(pp2)
}


#' @title A useful plot with NES and pvalue to display the result of GSEA analysis
#'
#' @param gseRes GSEA result generated by clusterProfiler.
#' @param UpperLim Display the 1 to <UpperLim> terms.
#' @param LowerLim Display the <LowerLim> to the final result. It is usually the total number minus the required number and then adding one.
#' @param Col1 Color 1
#' @param Col2 Color 2
#' @param SepCol Method of display result. Display according groups or P values. if True display by groups, false displayed by pval.default is false.
#' @param HiExprGroupName Char vector, when SepCol=T, define the hi-expr group name.
#' @param LoExprGroupName Char vector, when SepCol=T, define the low-expr group name.
#' @param WarpWidth Wrap length
#'
#' @return a ggplot2 object.
#' @export
#' @author Zhiming Ye
#'
#' @examples
PlotNES<-function(gseRes,UpperLim,LowerLim,Col1="#B39CD0",Col2="#FBEAFF",SepCol=F,HiExprGroupName="A",LoExprGroupName="B",WarpWidth=50){
  plotGSEA<-as.data.frame(gseRes)%>%dplyr::arrange(desc(NES))
  cat("There are ",length(plotGSEA$NES)," items. ",sum(plotGSEA$NES>0)," terms may upregulate, and ",sum(plotGSEA$NES<0)," terms may downregulate.\n")
  plotGSEA<-plotGSEA[c(c(1:UpperLim),c(LowerLim:nrow(plotGSEA))),]
  if(!SepCol)
  {
    p<-ggplot(plotGSEA, aes(NES, fct_reorder(Description, NES), fill=qvalues)) +geom_col(orientation='y') +scale_fill_continuous(low=Col1, high=Col2, guide=guide_colorbar(reverse=F)) + scale_y_discrete(labels=function(x) str_wrap(x, width=WarpWidth))+theme_bw() + ylab(NULL)
    return(p)}
  else{
    plotGSEA%<>%dplyr::mutate(Enrichment=ifelse(NES>0,HiExprGroupName,LoExprGroupName))
    plotGSEA$Enrichment<-factor(plotGSEA$Enrichment,levels = c(HiExprGroupName,LoExprGroupName))
    p<-ggplot(plotGSEA, aes(NES, fct_reorder(Description, NES), fill=Enrichment)) +geom_col(orientation='y') +scale_fill_manual(values=c(Col1,Col2)) + scale_y_discrete(labels=function(x) str_wrap(x, width=WarpWidth))+theme_bw() + ylab(NULL)
    return(p)
  }
}

#' @title Remove "_" in MsigDB items
#'
#' @param CHAR vector of MsigDB item names
#' @param start number of character which is used to start to separate.
#'
#' @return
#' @export
#' @author Zhiming Ye
#'
#' @examples
PrettyMsigdbName<-function(CHAR,start){
  #library(Hmisc)
  #library(tidyverse)
  CHAR2<-capitalize(tolower(
    gsub(
      pattern = "_",
      replacement = " ",
      substr(CHAR, start, nchar(CHAR))
    )))
  return(CHAR2)
}
#' @title Generate Sankey Diagram.
#'
#' @param InfoTable A table with specific 2 columns which brings Categorical variable.
#' @param ItemA Column name of item A.
#' @param ItemB Column name of item B.
#' @param ItemAlist specified ordering, a character vector.
#' @param ItemBlist specified ordering, a character vector.
#'
#' @return A ggplot2 object
#' @export
#' @author Zhiming Ye
#'
#' @examples
PlotSanky <- function(InfoTable,
                      ItemA,
                      ItemB,
                      ItemAlist,
                      ItemBlist) {
  #library(tidyverse)
  #library(ggalluvial)
  InfoTable <-
    InfoTable %>% dplyr::select(!!sym(ItemA),!!sym(ItemB)) %>% na.omit()
  SankyInfo <-
    InfoTable %>% dplyr::mutate(Sanky = paste0(!!sym(ItemA), "《",!!sym(ItemB)))
  SankyTable <-
    data.frame(Items = names(table(SankyInfo$Sanky)), counts = as.numeric(table(SankyInfo$Sanky)))
  SankyTable <-
    SankyTable %>% dplyr::mutate(
      ItemsA = sapply(strsplit(Items, '《'), function(x)
        x[1]),
      ItemsB = sapply(strsplit(Items, '《'), function(x)
        x[2])
    )
  SankyTable$ItemsA <- factor(SankyTable$ItemsA, levels = ItemAlist)
  SankyTable$ItemsB <- factor(SankyTable$ItemsB, levels = ItemBlist)
  p <- ggplot(SankyTable,
              aes(y = counts,
                  axis1 = ItemsA, axis2 = ItemsB)) +
    geom_alluvium(
      aes(fill = ItemsA),
      width = 0,
      reverse = FALSE,
      discern = TRUE
    ) + #控制线条流向
    geom_stratum(width = 0.07,
                 reverse = FALSE,
                 discern = F) + #控制中间框的宽度
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      reverse = FALSE,
      size = 4,
      angle = 0,
      discern = TRUE
    ) + #定义中间的文字
    scale_x_continuous(breaks = 1:2, labels = c(ItemA, ItemB)) + #定义X轴上图标排序
    theme_bw() +
    theme(legend.position = "none")
  return(p)
}
#' @title Generate pie plot
#'
#' @param dt A data frame contain items and frequency
#' @param itemAndfreq Vector specify the items and frequency, should be numeric.
#' @param pretty
#' @param Tag The label
#' @param Order Order, a character vector
#' @param col_type npg or lancet or nejm
#'
#' @return A ggplot2 object
#' @export
#' @author Zhiming Ye
#'
#' @examples
PlotPie <-
  function(dt,
           itemAndfreq = c(1, 2),
           pretty = T,
           Tag = "B",
           Order = "no",
           col_type = "nejm") {
    dt <- as.data.frame(dt)
    colnames(dt)[itemAndfreq[2]] <- "A"
    colnames(dt)[itemAndfreq[1]] <- "B"
    #library(ggpubr)
    #library(ggsci)
    dt = dt[order(dt$A, decreasing = TRUE), ]
    myLabel = as.vector(dt$B)
    myLabel = paste(myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)", sep = "")
    colnames(dt)[itemAndfreq[1]] <- "Tag"
    if (Order != "no") {
      dt$Tag <- factor(dt$Tag, levels = Order)
    }
    colnames(dt)[itemAndfreq[1]] <- Tag
    if (pretty) {
      if (col_type == "npg") {
        print(ggpie(dt, "A", fill = Tag, label = myLabel) + scale_fill_npg())
      }
      if (col_type == "lancet") {
        print(ggpie(dt, "A", fill = Tag, label = myLabel) + scale_fill_lancet())
      }
      if (col_type == "nejm") {
        print(ggpie(dt, "A", fill = Tag, label = myLabel) + scale_fill_nejm())
      }
    }
    else{
      ggpie(dt, "A", fill = Tag, label = myLabel)
    }
  }
#' @title Arrange Table Function - Sort and filter columns of a data frame or matrix by a specified vector
#'
#' @param Mat Target matrix or df.
#' @param FilterList Accord which to sort and filter columns
#'
#' @return a Matrix
#' @export
#' @author Zhiming Ye
#'
#' @examples
Arrange_Table <- function(Mat, FilterList) {
  cat("Should have rowname.\n")
  Mat <-
    as.data.frame(Mat) %>% dplyr::filter(rownames(Mat) %in% FilterList)
  rn1 <- rownames(Mat)
  rn2 <- FilterList
  Mat <- Mat[rn1[match(rn2, rn1)], ]
  return(Mat)
}
#' @title Do ssGSEA analysis
#'
#' @param Mat Expression Matrix
#' @param GMTfilePath gmt file path from MsigDB and so on.
#' @param numOfCore use how many CPU cores to calculate.
#'
#' @return A enrichment score matrix
#' @export
#' @author Zhiming Ye
#'
#' @examples
do.ssGSEA <- function(Mat, GMTfilePath,numOfCore=2) {
  parallel.sz = numOfCore
  parallel.type = 'SOCK'
  #library(GSVA)
  #library(GSEABase)
  genesetHypo <- getGmt(file.path(GMTfilePath))
  CGGA_gsvaHypo <-
    gsva(
      as.matrix(Mat),
      genesetHypo,
      method = 'ssgsea',
      kcdf = 'Gaussian',
      abs.ranking = TRUE,
      min.sz = 1,
      parallel.sz = parallel.sz
    )
  return(CGGA_gsvaHypo)
}
#' @title Do GSVA analysis
#'
#' @param Mat Expression Matrix
#' @param GMTfilePath gmt file path from MsigDB and so on.
#'
#' @return A GSVA score matrix
#' @export
#' @author Zhiming Ye
#'
#' @examples
do.GSVA <- function(Mat, GMTfilePath) {
  #library(GSVA)
  #library(GSEABase)
  genesetHypo <- getGmt(file.path(GMTfilePath))
  CGGA_gsvaHypo <-
    gsva(
      as.matrix(Mat),
      genesetHypo,
      method = 'gsva',
      kcdf = 'Gaussian',
      abs.ranking = TRUE,
      min.sz = 1
    )
  return(CGGA_gsvaHypo)
}
#' @title remove the Ensembl ID
#'
#' @param Mat Expression Matrix from Xena and so on
#' @param ColNum Which column contains ensembl ID, as numeric.
#'
#' @return A data frame
#' @export
#' @author Zhiming Ye
#'
#' @examples
rmEnsemblDot <- function(Mat, ColNum) {
  cat("New Col: SimSymbol.\n")
  Mat <-
    Mat %>% dplyr::mutate(SimSymbol = sapply(strsplit(Mat[, ColNum], '[.]'), function(x)
      x[1]))
  return(Mat)
}
#' @title remove separate from character
#'
#' @param Chr Target character vector
#' @param SepChr seperate character
#' @param num Target character location
#'
#' @return a character vector
#' @export
#' @author Zhiming Ye
#'
#' @examples
rmSep <- function(Chr,SepChr,num) {
  SimSymbol<-sapply(strsplit(Chr,SepChr), function(x)x[num])
  return(SimSymbol)
}
#' @title Convert Ensembl ID to Symbol
#'
#' @param TargetGene Target gene ID
#' @param ColNum Number of column contains ensembl ID, only nesessary when input a matrix
#' @param OrgDB The annotation DB
#' @param Type to:ensembl ID->symbol, from:symbol->ID
#'
#' @return a matrix
#' @export
#' @author Zhiming Ye
#'
#' @examples
MapSymbol <- function(TargetGene, ColNum=1, OrgDB,Type="to") {
  if(typeof(TargetGene)=="character"){
    if(Type=="to"){
      mapIds(x = OrgDB, keys = TargetGene,column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    }else{
      mapIds(x = OrgDB, keys = TargetGene,column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
    }

  }else{
    if(Type=="to"){
      cat("New Col: Symbol. \n")
      #library(org.Hs.eg.db)
      Mat$Symbol <- mapIds(
        x = org.Hs.eg.db,
        keys = Mat[, ColNum],
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
      return(Mat)
    }else{
      cat("New Col: Symbol. \n")
      #library(org.Hs.eg.db)
      Mat$Symbol <- mapIds(
        x = org.Hs.eg.db,
        keys = Mat[, ColNum],
        column = "ENSEMBL",
        keytype = "SYMBOL",
        multiVals = "first"
      )
      return(Mat)
    }
  }

}

#' @title NMF Estimation
#'
#' @param mat Matrix
#' @param pdfname Output PDF location
#' @param pngname Output PNG location
#' @param filename Output RData location
#'
#' @return NULL
#' @export
#' @author Zhiming Ye
#'
#' @examples
NMFEstimate <- function(mat, pdfname, pngname, filename) {
  #library(NMF)
  estimate <- nmf(
    mat,
    rank = 2:10,
    method = "brunet",
    nrun = 50,
    seed = 486,
    .options="p6"
  )
  save(estimate, file = filename)
  pdf(pdfname, width = 8, height = 6)
  try(print(plot(estimate)))
  dev.off()
  png(pngname,
      width = 4400,
      height = 4000,
      res = 300)
  try(consensusmap(
    estimate,
    annRow = NA,
    annCol = NA,
    main = "Consensus matrix",
    info = FALSE
  ))
  dev.off()
}

#' @title Calculate the SFT value for WGCNA analysis
#'
#' @param Mat Expression matrix
#' @param OutPutFile Output figure location
#'
#' @return NULL
#' @export
#' @author Zhiming Ye
#'
#' @examples
WGCNACalcSFT <- function(Mat, OutPutFile) {
  #library(WGCNA)
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft = pickSoftThreshold(Mat, powerVector = powers, verbose = 5)
  cat(paste0("sft PowerEstimate is:", sft$powerEstimate, "\n"))
  PowerEst <- sft$powerEstimate
  pdf(OutPutFile, width = 8, height = 5)
  par(mfrow = c(1, 2))
  cex1 = 0.9
  plot(
    sft$fitIndices[, 1],-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence")
  )

  text(
    sft$fitIndices[, 1],-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )

  abline(h = 0.90, col = "red")
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
  )

  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  dev.off()
  return(PowerEst)
}


#' @title Output expression matrix based on top MAD value
#'
#' @param mat Expression matrix
#' @param num select how many genes
#'
#' @return A matrix
#' @export
#' @author Zhiming Ye
#'
#' @examples
CalcMad <- function(mat, num) {
  mads <- apply(mat, 1, mad)
  mat <- mat[rev(order(mads)), ]
  return(mat[1:num, ])
}

#' @title Batch uni-Cox regression to find survival-related genes
#'
#' @param StatusAndTimeCol Numeric vector, item first is number of status column and the second is the time column
#' @param DATA Data set
#'
#' @return A Cox-regression result
#' @export
#' @author Zhiming Ye
#'
#' @examples
UniCox <- function(StatusAndTimeCol, DATA) {
  # cat("STATUS=OS,TIME=OStime. Argument 1 is a collection.\n")
  covariates <- colnames(DATA)[-StatusAndTimeCol]
  univ_formulas <- sapply(covariates,
                          function(x)
                            as.formula(paste('Surv(OStime, OS)~', x)))
  #library(survival)
  #library(survminer)
  univ_models <-
    lapply(univ_formulas, function(x) {
      coxph(x, data = DATA)
    })
  univ_results <- lapply(univ_models,
                         function(x) {
                           x <- summary(x)
                           #获取p值
                           p.value <-
                             signif(x$wald["pvalue"], digits = 5)
                           #获取HR
                           HR <- signif(x$coef[2], digits = 5)

                           #获取95%置信区间
                           HR.confint.lower <-
                             signif(x$conf.int[, "lower .95"], 5)
                           HR.confint.upper <-
                             signif(x$conf.int[, "upper .95"], 5)
                           res <- c(p.value, HR,HR.confint.lower,HR.confint.upper)
                           names(res) <-
                             c("p.value", "HR","L","H")
                           return(res)
                         })
  #转换成数据框，并转置
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  return(res)
}

#' @title Perform a normality test by batch
#'
#' @param Table Target table
#'
#' @return P Val
#' @export
#' @author Zhiming Ye
#'
#' @examples
NormalizeTest <- function(Table) {
  MitoAndHypo2 <- as.matrix(Table)
  for (i in 1:ncol(MitoAndHypo2)) {
    print(shapiro.test(MitoAndHypo2[, i]))
  }

}

#' @title Write tab table file
#'
#' @param table A character, name of data frame or matrix object.
#' @param withRowname with row names or not.
#'
#' @return NULL, save a file to work dictionary.
#' @export
#' @author Zhiming Ye
#'
#' @examples
Write.tables <- function(table, withRowname = T) {
  table <- as.character(table)
  # cat("will save to the working dictionary.\n")
  PCAmatnolog <- get(table)
  if (withRowname) {
    PCAmatnolog <-
      PCAmatnolog %>% as.data.frame() %>% rownames_to_column(var = "temp3wr4ts")
    colnames(PCAmatnolog)[1] <- ""
  }
  write.table(
    PCAmatnolog,
    file = paste0(table, ".txt"),
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
  )
}

#' @title Convert FPKM to TPM
#'
#' @param expMatrix Expression matrix
#'
#' @return TPM matrix
#' @export
#' @author Zhiming Ye
#'
#' @examples
fpkmToTpm <- function(expMatrix)
{
  cat("Might be not correct\n")
  fpkmToTpm2 <- function(fpkm)
  {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  tpms <- apply(expMatrix, 2, fpkmToTpm2)
  return(tpms)
}

#' @title Estimate correlation by batch
#'
#' @param genegsvamatrix Target matrix or data frame.
#' @param TargetSet Target estimated column name
#' @param UseMethod default as "spearman"
#'
#' @return R and correlation coeffection
#' @export
#' @author Zhiming Ye
#'
#' @examples
CorrEstimate0 <-
  function(genegsvamatrix,
           TargetSet = "HALLMARK_HYPOXIA",
           UseMethod = "spearman") {
    data <- genegsvamatrix
    gene_name1 <- c()
    gene_name2 <- c()
    cor_r <- c()
    pvalue <- c()
    i = which(rownames(data) == TargetSet)
    for (r in 1:nrow(data)) {
      g1 = rownames(data)[i]
      g2 = rownames(data)[r]
      c_r = cor(as.numeric(data[i, ]), as.numeric(data[r, ]), method = UseMethod)
      p = cor.test(as.numeric(data[i, ]), as.numeric(data[r, ]), method =
                     UseMethod)[[3]]
      gene_name1 = c(gene_name1, g1)
      gene_name2 = c(gene_name2, g2)
      cor_r = c(cor_r, c_r)
      pvalue = c(pvalue, p)
      cat(paste0(r, "   done!\n"))

    }
    data_cor <- data.frame(gene_name1, gene_name2, cor_r, pvalue)
    return(data_cor)
  }

#' @title Build up a DEG data frame for clusterProfiler multi-group ORA analysis
#' @description Convert symbol, and add up group name in column.
#' @param Gene A vector containing genes.
#' @param Group A character of group name
#' @param logFC The vector of fold-change value
#' @param OrgDB Default as org.Hs.eg.db
#'
#' @return A data frame. Following rbind command can help you build up the full DEG list.
#' @export
#' @author Zhiming Ye
#'
#' @examples
BuildMultigroupDEGlist<-function(Gene,Group,logFC,OrgDB=org.Hs.eg.db){
  Genetable <- data.frame(Gene = Gene, Group = Group,Type=ifelse(logFC>0,"Up","Down"))
  ENTREZIDtable <- clusterProfiler::bitr(Genetable$Gene, fromType = "SYMBOL",
                                         toType = "ENTREZID", OrgDb = OrgDB)
  colnames(ENTREZIDtable)[1] <- "Gene"
  Genetable <- Genetable %>% left_join(ENTREZIDtable)
  return(Genetable)
}

.onAttach<-function(libname,pkgname){
  packageStartupMessage("\n***BulkSeqAssistant*** V5.5\n=============\nPart of the code came from online materials,\nFor learning and communication purposes only.\n=============\nAuthor:Zhiming Ye, SKLRD @ Guangzhou Medical University\n")
}


summary.PMCMR <- function(object, ...)
{
  OK <- inherits(object, c("PMCMR"))
  if (!OK)
    stop ("Not an object of class PMCMR")
  if (!is.matrix(object$statistic))
    stop ("Matrix object$statistic not found.")
  pval <- as.numeric(object$p.value)
  stat <- as.numeric(object$statistic)
  grp1 <- as.numeric(c(col(object$p.value)))
  cnam <- colnames(object$p.value)
  grp2 <- as.numeric(c(row(object$p.value)))
  rnam <- rownames(object$p.value)
  STAT <- object$dist

  if (!is.null(object$alternative)) {
    if (object$alternative == "less"){
      H0 <- paste(rnam[grp2], "-", cnam[grp1], ">=", "0")
      PVAL <- paste("Pr(<", STAT, ")", sep="")
    } else if (object$alternative == "greater"){
      H0 <- paste(rnam[grp2], "-", cnam[grp1], "<=", "0")
      PVAL <- paste("Pr(>", STAT, ")", sep="")
    } else {
      H0 <- paste(rnam[grp2], "-", cnam[grp1], "==", "0")
      PVAL <- paste("Pr(>|", STAT, "|)", sep="")
    }
  } else {
    H0 <- paste(rnam[grp2], "-", cnam[grp1], "==", "0")
    PVAL <- paste("Pr(>|", STAT, "|)", sep="")
  }

  STAT2 <- paste0(STAT, " value")
  OK <- !is.na(pval)
  ## Symbols
  symp <- symnum(pval[OK], corr=FALSE,
                 cutpoints = c(0,  .001,.01,.05, .1, 1),
                 symbols = c("***","**","*","."," "))

  xdf <- data.frame(statistic = round(stat[OK], 3),
                    p.value = format.pval(pval[OK]),
                    symp)
  rownames(xdf) <- H0[OK]
  names(xdf) <- c(STAT2, PVAL, "")
  ##
  return(xdf)
}


#' @title Multi group Parametric test
#'
#' @param Mat Columns should be every variable to be tested, and another one column contains information about grouping, which is refer to the FactorCol
#' @param FactorCol Which column contains the grouping information. The number of it.
#' @param Method One of "LSD" or "SNK"
#' @param OnlySig Whether show only significant result.
#' @param pcutoff P value cut off.
#'
#' @return
#' @export
#' @author Zhiming Ye
#'
#' @examples
MultiAOV<-function(Mat,FactorCol,Method="LSD",OnlySig=F,pcutoff=0.01){
  # library(PMCMRplus)
  WillTest<-Mat
  ResultRes<-data.frame()
  colnames(WillTest)[FactorCol]<-"Cluster"
  CellTypeList<-colnames(WillTest)[-FactorCol]
  for(CellType in CellTypeList){
    TestSet<-WillTest[,c("Cluster",CellType)]
    colnames(TestSet)[2]<-"Cell"
    # KwRes<-kruskal.test(Cell~Cluster,TestSet)
    TestSet$Cluster<-as.factor(TestSet$Cluster)
    QA<-aov(Cell~Cluster,TestSet)
    summary(QA)[[1]][["F value"]][1]->Fvalue
    if(Fvalue<pcutoff){
      if(Method!="LSD"){
        cat(paste0("ANOVA PASS : ",CellType," , Processing SNK Test...\n"))
        WRes<-snkTest(QA)
      }else{
        cat(paste0("ANOVA PASS : ",CellType," , Processing LSD Test...\n"))
        WRes<-lsdTest(QA)
      }
      Result0<-summary.PMCMR(WRes)[,-3]%>%as.data.frame()%>%rownames_to_column(var="Comparasion")%>%dplyr::mutate(GroupName=CellType)
      ResultRes<-rbind(ResultRes,Result0)
    }
  }
  if(OnlySig){
    ResultRes2<-ResultRes%>%dplyr::filter(P.adj<pcutoff)
    ResultRes3<-ResultRes%>%dplyr::filter(GroupName%in%names(table(ResultRes2$GroupName)))
    return(ResultRes3)
  }
  else{
    return(ResultRes)
  }
}


#' @title Extract Marker genes of scanpy DEG Result
#' @description The result is generated by pandas using for command. Ref to the scanpy document.
#' @description pd.DataFrame(
#' @description     {group + '_' + key[:1]: result[key][group]
#' @description     for group in groups for key in ['names', 'pvals']}).head(5)
#' @param DEGmat Result of DEG
#' @param QueryGene Gene Symbol
#' @param numOfClusters Number of cell clusters
#' @param verbose Print Pval and log2FC
#' @param PvalCutoff Pval cut off
#' @param lFCcutoff log2FC cut off
#' @param IsHighExpGenes Choose whether FC>0 or <0
#' @param geneColName
#' @param pvalColName
#' @param lfcColName
#'
#' @return Print the chosen cell cluster
#' @export
#' @author Zhiming Ye
#'
#' @examples
scanpy.QueryDEG<-function(DEGmat,QueryGene,numOfClusters,verbose=T,PvalCutoff=0.01,lFCcutoff=1,IsHighExpGenes=T,geneColName="_names",pvalColName="_pvals_adj",lfcColName="_logfoldchanges"){
  cat("*Query Gene : ",QueryGene,"\n")
  for(i in 0:numOfClusters){
    gene_c<-paste0(i,geneColName)
    pval_c<-paste0(i,pvalColName)
    lFC_c<-paste0(i,lfcColName)
    QueryGene<-paste0("^",QueryGene,"$")
    Location<-which(grepl(QueryGene,as.vector(unlist(DEGmat[,gene_c]))))
    if(IsHighExpGenes){
      if(unlist(DEGmat[Location,pval_c])<PvalCutoff&unlist(DEGmat[Location,lFC_c])>lFCcutoff){
        if(verbose){
          cat("Cell Cluster ",i,"     Pval: ",unlist(DEGmat[Location,pval_c])," , lFC: ",unlist(DEGmat[Location,lFC_c]),"\n")
        }else{
          cat("Cell Cluster ",i,"\n")
        }
      }
    }else{
      if(unlist(DEGmat[Location,pval_c])<PvalCutoff&unlist(DEGmat[Location,lFC_c])<lFCcutoff){
        if(verbose){
          cat("Cell Cluster ",i,"     Pval: ",unlist(DEGmat[Location,pval_c])," , lFC: ",unlist(DEGmat[Location,lFC_c]),"\n")
        }else{
          cat("Cell Cluster ",i,"\n")
        }
      }
    }
  }
}


#
#  ggbiplot.r
#
#  Copyright 2011 Vincent Q. Vu.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#
ggbiplot.internal <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                     obs.scale = 1 - scale, var.scale = scale,
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68,
                     labels = NULL, labels.size = 3, alpha = 1,
                     var.axes = TRUE,
                     circle = FALSE, circle.prob = 0.69,
                     varname.size = 3, varname.adjust = 1.5,
                     varname.abbrev = FALSE, ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)

  stopifnot(length(choices) == 2)

  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }

  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))

  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])

  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)

  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }

  # Scale the radius of the correlation circle so that it corresponds to
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)

  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))

  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }

  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs,
                       sprintf('(%0.1f%% explained var.)',
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))

  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }

  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }

  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }

  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)

  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) +
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()

  if(var.axes) {
    # Draw circle
    if(circle)
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'),
                         size = 1/2, alpha = 1/3)
    }

    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')),
                   color = muted('red'))
  }

  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups),
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)
    }
  }

  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))

    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }

  # Label the variable axes
  if(var.axes) {
    g <- g +
      geom_text(data = df.v,
                aes(label = varname, x = xvar, y = yvar,
                    angle = angle, hjust = hjust),
                color = 'darkred', size = varname.size)
  }
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)),
  #                               palette = 'Dark2')
  # }

  # TODO: Add a second set of axes

  return(g)
}
