#library(psych,warn.conflicts=F)
#library(stats,warn.conflicts=F)
#' Calculates measure M
#' 
#' This function calculates measure M according to algorithm of geNorm
#'
#' @param expression a matrix of expression levels. Each row corresponds to a sample and each column to a gene.
#' @param ctVal a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.
#' @return A sorted dataframe with two columns, 'Genes' and 'M' (the relative stability; lower means more stable).
#' @export
#' @references 
#' Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
#' @examples
#' FIB
#' measureM(FIB,ctVal=FALSE)
#' FIBct
#' measureM(FIBct)
#' @importFrom "stats" "sd"
measureM=function(expression,ctVal=TRUE)
{
  if(ctVal)
  {
    dct<-apply(expression, 2, function(x) {x-min(x)})
    expression=2^-dct
  }
  else
  {
    expression<-apply(expression, 2, function(x) {x/max(x)})
  }
  m=nrow(expression)
  n=ncol(expression)
  M=list()
  for(j in 1:n)
  {
    Vjk=list()
    for(k in 1:n)
    {
      Ajk=list()
      #if (j==k)next()
      for(i in 1:m)
      {
        Ajk=c(Ajk,log2(expression[i,j]/expression[i,k]))
      }
      Vjk=c(Vjk,sd(as.matrix(Ajk)))
      
    }
    M=c(M,sum(unlist(Vjk))/(n-1))
  }
  M_a<-data.frame(Genes=colnames(expression),M=unlist(M))
  return (M_a[order(M_a[,2],decreasing=T),])
}

#' Ranks genes
#' 
#' Uses the geNorm algorithm to determine the most stably expressed genes.
#'
#' @param expression a matrix of expression levels. Each row corresponds to a sample and each column to a gene.
#' @param genes a data frame to output the result of the function
#' @param ctVal a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.
#' @return A sorted dataframe with two columns, 'Genes' and 'Avg.M'. The last two genes are the two most stable control genes.
#' @return Avg.M is average expression stability values (M) of remaining control genes during stepwise exclusion of the least stable control gene.
#' @export
#' @references 
#' Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
#' @examples
#' FIB
#' geNorm(FIB,ctVal=FALSE)
#' FIBct
#' geNorm(FIBct)
geNorm=function(expression,genes=data.frame(Genes=character(0),Avg.M=numeric(0)),ctVal=TRUE)
{
  geNorm_result=measureM(expression,ctVal)
  #print(mean(geNorm_result$M))
  n=ncol(expression)
  if(n<=2)
  {
    stable_genes=paste(geNorm_result[1,1],geNorm_result[2,1], sep = "-")
    genes=rbind(genes,  data.frame(Genes=stable_genes,Avg.M=mean(geNorm_result$M)) )
    return (genes)
  }
  else
  {
    index<-which(colnames(expression)==geNorm_result[1,1])
    genes=rbind(genes,  data.frame(Genes=geNorm_result[1,1],Avg.M=mean(geNorm_result$M)) )
    geNorm(expression[,-index],genes,ctVal)
  }
}

#' Ranks genes
#' 
#' Uses the geNorm algorithm to determine the most stably expressed genes.
#'
#' @param expression a matrix of expression levels. Each row corresponds to a sample and each column to a gene.
#' @param genes a data frame to output the result of the function
#' @param ctVal a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.
#' @return A sorted dataframe with two columns, 'Genes' and 'Avg.M'. The last two genes are the two most stable control genes. 
#' @return Avg.M is average expression stability values (M) of remaining control genes during stepwise exclusion of the least stable control gene.
#' @export
#' @references 
#' Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
#' @examples
#' FIB
#' geNorm2(FIB,ctVal=FALSE)
#' FIBct
#' geNorm2(FIBct)
geNorm2=function(expression,genes=data.frame(Genes=character(0),Avg.M=numeric(0)),ctVal=TRUE)
{
  geNorm_result=measureM(expression,ctVal)
  #print(mean(geNorm_result$M))
  n=ncol(expression)
  if(n==2)
  {
    genes=rbind(genes,  data.frame(Genes=geNorm_result[1,1],Avg.M=mean(geNorm_result$M)) )
    genes=rbind(genes,  data.frame(Genes=geNorm_result[2,1],Avg.M=NA) )
    return (genes)
  }
  else
  {
    index<-which(colnames(expression)==geNorm_result[1,1])
    genes=rbind(genes,  data.frame(Genes=geNorm_result[1,1],Avg.M=mean(geNorm_result$M)) )
    geNorm2(expression[,-index],genes,ctVal)
  }
}

#' Calculates V(n+1/n) values
#' 
#' Useful for establishing the quality of your normalization regime. See Vandesompele 2002 for advice on interpretation.
#'
#' @param expression a matrix of expression levels. Each row corresponds to a sample and each column to a gene.
#' @param ctVal a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.
#' @return A Series of values [V2/3, V3/V4, V4/V5, ...]. 
#' @export
#' @references 
#' Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
#' @examples
#' FIB
#' pairwiseV(FIB,ctVal=FALSE)
#' FIBct
#' pairwiseV(FIBct)
#' @importFrom  "stats" "sd"
pairwiseV=function(expression,ctVal=TRUE)
{
  buf=geNorm(expression,ctVal=F)
  buf=buf$Genes
  genes=strsplit(as.character(buf[length(buf)]),split = "-")
  n=length(buf)
  for(i in 1:(n-1))
  {
    genes=c(genes,as.character(buf[n-i]))
  }
  genes=unlist(genes)
  expression=as.data.frame( expression[,genes])
  
  Vs=data.frame(V=character(0),Value=numeric(0))
  for(i in 2:n) 
  {
    genesn=as.matrix(expression[,c(1:i)])
    NFn=apply(genesn, 1, psych::geometric.mean)
    
    genesn1=as.matrix(expression[,c(1:(i+1))])
    NFn1=apply(genesn1, 1, psych::geometric.mean)
    
    An_n1=log2(NFn/NFn1)
    Vn_n1=sd(An_n1)
    Vs=rbind(Vs, data.frame(V=paste("V",i,"/","V",i+1,sep = ""),Value=Vn_n1) )
  }
  
  return (Vs)
}

#' Plots average M of remaining genes
#' 
#' This function plots the average expression stability values of remaining control genes.
#'
#' @param Mrem the result returned by function of geNorm()
#' @return NULL
#' @export
#' @references 
#' Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
#' @examples
#' FIB
#' x=geNorm(FIB,ctVal=FALSE)
#' plotM(x)
#' 
#' FIBct
#' y=geNorm(FIBct)
#' plotM(y)
#' @importFrom "graphics" "axis" "barplot" "par" "plot" "text" "title" 
plotM=function(Mrem)
{
  par(mar = rep(5.1,4.1,4.1,2.1),mgp=c(4,1,0))
  x<-c(0.5:(nrow(Mrem)-0.5))
  plot(x,Mrem$Avg.M ,type="o",col="blue",frame.plot=TRUE,xaxt='n',pch=16,ylim = c(0,max(Mrem$Avg.M)),lwd=1,ylab="",xlab="")
  xlabs=Mrem$Genes
  lastgenes=gsub(pattern = "-", replacement = '\n', x = as.character(xlabs[length(xlabs)]))   
  xlabs=c(as.character(xlabs[-length(xlabs)]),lastgenes)
  axis(1,x,xlabs,las=2)
  Title='Average expression stability values of remaining control genes'
  yLabel='Average expression stability M'
  xLabel='<:::::  Least stable genes             Most stable genes ::::>'
  
  title(main=Title, ylab=yLabel, xlab=xLabel) 
}

#' Plots V(n+1/n) values
#' 
#' This function plots the average expression stability values of remaining control genes.
#'
#' @param Vs the result returned by function of pairwiseV()
#' @return NULL
#' @export
#' @references 
#' Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
#' @examples
#' FIB
#' Vs1=pairwiseV(FIB,ctVal=F)
#' plotV(Vs1)
#' 
#' FIBct
#' Vs2=pairwiseV(FIBct)
#' plotV(Vs2)
#' @importFrom "graphics" "par" "barplot" "text" 
plotV=function(Vs)
{
  par(mar = rep(5.1,4.1,4.1,2.1),mgp=c(3,1,0))
  Title='Determination of the optimal number of control genes for normalization'
  barplot(Vs$Value,names.arg= Vs$V,xlab="Pairwise Variations",main=Title,ylim=c(0,max(Vs$Value)),col = "yellow")
  text(1:(nrow(Vs))*1.2-.5,Vs$Value/2, round(Vs$Value,3))
}

#' The normalized expression level of the ten housekeeping genes in fibroblast cells
#' @name FIB
#' @aliases FIB
#' @title Reload Saved Datasets
#' @description The normalized expression level of the ten housekeeping genes in fibroblast cells
#' @references 
#' Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
NULL

#' The CT values of the ten housekeeping genes in fibroblast cells
#' @name FIBct
#' @aliases FIBct
#' @title Reload Saved Datasets
#' @description The CT values of the ten housekeeping genes in fibroblast cells
#' @references 
#' Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
NULL