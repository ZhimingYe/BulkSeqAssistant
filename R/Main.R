#Last changed : 2022-9-2 11:20
.onAttach<-function(libname,pkgname){
  packageStartupMessage("***BulkSeqAssistant*** V5.1.6\nAuthor:Zhiming Ye, Guangzhou Medical University\nPart of the code comes from online materials for learning and communication purposes only.\nBug feedback Email: zhiming.ye@qq.com\n")
}
#library(tidyverse)
#library(magrittr)
#' @title Get DEG from DEseq2 result
#'
#' @param DEseq2Result A DEseq2 result. Should be converted to data frame.
#' @param Type Can be "all", "up" and "down".
#' @param Padj P adjusted cut off value.
#' @param lFC logFoldChange cut off value.
#'
#' @return A filtered DEG data frame.
#' @export
#'
#' @examples
getDESeq2DEG<-function(DEseq2Result,Type,Padj=0.05,lFC=1){
  cat("Type can be all or up or down.\n")
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
#'
#' @examples
wilDEGs <- function(expr_matrix, sample_df, pval, fdr, FC, adj.method = "BH"){
  # tot_sam_ct <- ncol(expr_matrix)
  # tot_gene_ct <- nrow(expr_matrix)
  # tot_sam_label <- nrow(sample_df)
  cat("input1:data_matrix，rowname is gene，colname is sample\ninput2:sample_group(df),column must contains Sample_ID and Type\n")
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
#' @param Gene
#' @param log2FC
#' @param FromType default is "SYMBOL"
#'
#' @return A ranked named numeric vector. Names of the numbers is the ENTREZID.
#' @export
#'
#' @examples
Ranked.GS<-function(Gene,log2FC,FromType = "SYMBOL"){
  Genetable<-data.frame(Gene=Gene,log2FC=log2FC)
  #library(clusterProfiler)
  #library(enrichplot)
  #library(ReactomePA)
  #library(org.Hs.eg.db)
  ENTREZIDtable<-clusterProfiler::bitr(Genetable$Gene,fromType = FromType,toType = "ENTREZID",OrgDb = org.Hs.eg.db)
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
#'
#' @return An enrichment result generated by clusterProfiler.
#' @export
#'
#' @examples
doGO<-function(GS,ont="BP",PVal=0.01,QVal=0.05){
  #library(clusterProfiler)
  #library(enrichplot)
  #library(ReactomePA)
  #library(org.Hs.eg.db)
  GeneList<-GS%>%unique()%>%bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  EnrichGO<-enrichGO(GeneList$ENTREZID,ont = ont,OrgDb = org.Hs.eg.db,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable(OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  return(EnrichGO)
}
#' @title ORA KEGG analysis by clusterProfiler by symbol identifiers.
#'
#' @param GS Gene Symbols for enrichment analysis. NOT ENSEMBL ID.
#' @param PVal P adjusted value cut off, default 0.01
#' @param QVal qvalue cut off, default 0.05
#'
#' @return An enrichment result generated by clusterProfiler.
#' @export
#'
#' @examples
doKEGG<-function(GS,PVal=0.01,QVal=0.05){
  #library(clusterProfiler)
  #library(enrichplot)
  #library(ReactomePA)
  #library(org.Hs.eg.db)
  GeneList<-GS%>%unique()%>%bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  EnrichKEGG<-enrichKEGG(GeneList$ENTREZID,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable(OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  return(EnrichKEGG)
}
#' @title ORA Reactome analysis by clusterProfiler by symbol identifiers.
#'
#' @param GS Gene Symbols for enrichment analysis. NOT ENSEMBL ID.
#' @param PVal P adjusted value cut off, default 0.01
#' @param QVal qvalue cut off, default 0.05
#'
#' @return An enrichment result generated by clusterProfiler.
#' @export
#'
#' @examples
doRA<-function(GS,PVal=0.01,QVal=0.05){
  #library(clusterProfiler)
  #library(enrichplot)
  #library(ReactomePA)
  #library(org.Hs.eg.db)
  GeneList<-GS%>%unique()%>%bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  EnrichRA<-enrichPathway(GeneList$ENTREZID,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable(OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
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
  pic1<-ggbiplot(res.pca2,var.axes = F,obs.scale = 0.5,groups = as.factor(group),ellipse = Ellipse,circle = F)+theme_bw()
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
#' @param gseRes GSEA result generated by clusterProfiler
#' @param showItemsNum Numbers of shown items.
#' @param Col1 Color 1
#' @param Col2 Color 2
#' @param Ratio Ration of numbers of up regulated pathways to down regulated pathways.
#' @param SepCol Method of display result. Display according groups or P values. if True display by groups, false displayed by pval.default is false.
#' @param HiExprGroupName Char vector
#' @param LoExprGroupName Char vector
#' @param WarpWidth Wrap length
#'
#' @return a ggplot2 object.
#' @export
#'
#' @examples
PlotNES<-function(gseRes,showItemsNum,Col1="#B39CD0",Col2="#FBEAFF",Ratio=0.5,SepCol=F,HiExprGroupName="A",LoExprGroupName="B",WarpWidth=50){
  if(Ratio>=1|(!is.numeric(Ratio))){
    cat("ERROR!\n")
  }
  else{
    plotGSEA<-as.data.frame(gseRes)%>%dplyr::arrange(desc(NES))
    plotGSEA<-plotGSEA[c(c(1:round(showItemsNum*Ratio)),c(nrow(plotGSEA)-round(showItemsNum*(1-Ratio))):nrow(plotGSEA)),]
    if(!SepCol)
    {
      p<-ggplot(plotGSEA, aes(NES, fct_reorder(Description, NES), fill=qvalues)) +geom_col(orientation='y') +scale_fill_continuous(low=Col1, high=Col2, guide=guide_colorbar(reverse=TRUE)) + scale_y_discrete(labels=function(x) str_wrap(x, width=WarpWidth))+theme_bw() + ylab(NULL)
      return(p)}
    else{
      plotGSEA%<>%dplyr::mutate(Enrichment=ifelse(NES>0,HiExprGroupName,LoExprGroupName))
      plotGSEA$Enrichment<-factor(plotGSEA$Enrichment,levels = c(HiExprGroupName,LoExprGroupName))
      p<-ggplot(plotGSEA, aes(NES, fct_reorder(Description, NES), fill=Enrichment)) +geom_col(orientation='y') +scale_fill_manual(values=c(Col1,Col2)) + scale_y_discrete(labels=function(x) str_wrap(x, width=WarpWidth))+theme_bw() + ylab(NULL)
      return(p)
    }
  }
}

#' @title Remove "_" in MsigDB items
#'
#' @param CHAR vector of MsigDB item names
#' @param start number of character which is used to start to separate.
#'
#' @return
#' @export
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
#'
#' @examples
rmSep <- function(Chr,SepChr,num) {
  SimSymbol<-sapply(strsplit(Chr,SepChr), function(x)x[num])
  return(SimSymbol)
}
#' @title Convert Ensembl ID to Symbol
#'
#' @param Mat Target expression matrix
#' @param ColNum Number of column contains ensembl ID
#'
#' @return a matrix
#' @export
#'
#' @examples
MapSymbol <- function(Mat, ColNum) {
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
#'
#' @examples
UniCox <- function(StatusAndTimeCol, DATA) {
  cat("STATUS=OS,TIME=OStime. Argument 1 is a collection.\n")
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
#'
#' @examples
Write.tables <- function(table, withRowname = T) {
  table <- as.character(table)
  cat("will save to the working dictionary.\n")
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
#'
#' @examples
CorrlatinEstimate0 <-
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

