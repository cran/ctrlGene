#library(psych,warn.conflicts=F)
#library(stats,warn.conflicts=F)

#' Calculates descriptive statistics
#' 
#' This function calculates descriptive statistics of genes.
#'
#' @param expression a matrix of expression levels. Each row corresponds to a sample and each column to a gene.
#' @param ctVal a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.
#' @return A matrix of descriptive statistics: 
#' @return N: number of samples;
#' @return GM[CP]: the geometric mean of CP; 
#' @return AM[CP]: the arithmetic mean of CP; 
#' @return Min[CP] and Max [CP]: the extreme values of CP; 
#' @return SD[+/- CP]: the standard deviation of the CP; 
#' @return CV[CP]: the coefficient of variance expressed as a percentage on the CP level; 
#' @return Min[x-fold] and Max [x-fold]: the extreme values of expression levels expressed as an absolute x-fold over- or under-regulation coefficient; 
#' @return SD[+/- x-fold]: standard deviation of the absolute regulation coefficients.
#' @export
#' @references 
#' Pfaffl MW, Tichopad A, Prgomet C, Neuvians TP. Biotechnol Lett (2004) <doi: 10.1023/B:BILE.0000019559.84305.47>
#' @examples
#' FIBct
#' cpSta(FIBct)
cpSta=function(expression,ctVal=TRUE)
{
  if (!ctVal){expression=log2(expression)}
  N=rep(nrow(expression), times=ncol(expression))
  GM_CP=apply(expression, 2, psych::geometric.mean)
  AM_CP=apply(expression, 2, mean)  
  Min_CP=apply(expression, 2, min) 
  Max_CP=apply(expression, 2, max) 
  
  AVEDEV=function(x)
  {
    return(sum(abs(x-mean(x)))/length(x))
  }
  SD_CP=apply(expression, 2, AVEDEV)
  CV_CP=100*SD_CP/AM_CP
  
  Min_x_fold=-1/2^(Min_CP-GM_CP)
  Max_x_fold=2^(Max_CP-GM_CP)
  
  SD_X_fold_f=function(x)
  {
    y=ifelse(x >mean(x),x-mean(x),mean(x)-x)
    return(AVEDEV(2^y))
  }
  SD_x_fold=apply(expression, 2, SD_X_fold_f) 
  
  rz=matrix(c(N,GM_CP,AM_CP,Min_CP,Max_CP,SD_CP,CV_CP,Min_x_fold,Max_x_fold,SD_x_fold),ncol=length(GM_CP),nrow = 10,byrow=T)
  colnames(rz)=colnames(expression)
  rownames(rz)=c("N","GM[CP]","AM[CP]","Min[CP]","Max[CP]","SD[+/- CP]","CV[%CP]","Min[x-fold]","Max[x-fold]","SD[+/- x-fold]")
  return(round(rz,3))
}


#' Analyzes pair-wise correlation
#' 
#' This function performs numerous pair-wise correlation analyses of genes. 
#' Within each such correlation the Pearson correlation coefficient (r) and the probability p value are calculated.
#' 
#' @param expression a matrix of expression levels. Each row corresponds to a sample and each column to a gene.
#' @param ctVal a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.
#' @return A matrix of the Pearson correlation coefficient (r) and the probability p value.
#' @export
#' @references 
#' Pfaffl MW, Tichopad A, Prgomet C, Neuvians TP. Biotechnol Lett (2004) <doi: 10.1023/B:BILE.0000019559.84305.47>
#' @examples
#' FIBct
#' pearsonCor(FIBct)
#' @importFrom "stats" "cor.test"
pearsonCor=function(expression,ctVal=TRUE)
{
  if (!ctVal){expression=log2(expression)}
  n=ncol(expression)
  rz=matrix(ncol=n,nrow=2*(n-1))
  colnames(rz)=colnames(expression)
  rowname_rz=list()
  for(i in 1:n)
  {
    for(j in 2:n) 
    {
      COT_TEST=cor.test(expression[,i],expression[,j])
      rz[2*(j-1)-1,i]=COT_TEST$estimate
      rz[2*(j-1),i]=COT_TEST$p.value
    }
    if (i+1>n) break;
    rowname_rz=c(rowname_rz,colnames(expression[i+1]),'p-value')
  }
  rownames(rz)=unlist(rowname_rz)
  return(round(rz,3))
}

#' Analyzes genes versus BestKeeper index
#' 
#' All genes are combined into an index. Then, correlation between each genes and the index is calculated, describing the relation between the index and the contributing gene by the Pearson correlation coefficient (r), coefficient of determination (r2) and the p-value.
#' 
#' @param expression a matrix of expression levels. Each row corresponds to a sample and each column to a gene.
#' @param ctVal a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.
#' @return A matrix of the Pearson correlation coefficient (r), coefficient of determination (r2) and the p-value.
#' @export
#' @references 
#' Pfaffl MW, Tichopad A, Prgomet C, Neuvians TP. Biotechnol Lett (2004) <doi: 10.1023/B:BILE.0000019559.84305.47>
#' @examples
#' FIBct
#' bki(FIBct)
#' @importFrom "stats" "cor.test"
bki=function(expression,ctVal=TRUE)  
{
  if (!ctVal){expression=log2(expression)}
  BKI=apply(expression, 1, psych::geometric.mean)
  n=ncol(expression)
  rz=matrix()
  for(i in 1:n)
  {
    COT_TEST=cor.test(expression[,i],BKI)
    r2=COT_TEST$estimate^2
    count=floor((length(expression[,i])+length(BKI))/2)
    
    SC_xy=sum(expression[,i]*BKI)-sum(expression[,i])*sum(BKI)/count  #HU113
    SS_x=sum(BKI^2)-sum(BKI)^2/length(BKI)#HS113
    SS_y= sum(expression[,i]^2)-sum(expression[,i])^2/length(expression[,i])#HT113
    slope=SC_xy/SS_x
    intercept=mean(expression[,i])-slope* mean(BKI)
    SE=(SS_y*(1-r2)/(count-2))^0.5
    Power_x_fold=2^slope
    
    buf=matrix(c(COT_TEST$estimate,r2,intercept,slope,SE,COT_TEST$p.value,Power_x_fold))
    
    if (length(rz)==1){rz=buf}  
    else {rz=cbind(rz,as.matrix( buf))}
  }
  colnames(rz)=colnames(expression)
  rownames(rz)=c("coeff. of corr. [r]","coeff. of det. [r^2]","intercept [CP]","slope [CP]","SE [CP]","p-value","Power [x-fold]")
  return(round(rz,3))
}

#' Determines stability of genes
#' 
#' This function combines the results of cpSta(), pearsonCor() and bki().
#' 
#' @param expression a matrix of expression levels. Each row corresponds to a sample and each column to a gene.
#' @param ctVal a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.
#' @return A list containing CP.statistics, pair.Wise.cor and HKG.vs.BestKeeper, which are returned by cpSta(), pearsonCor() and bki(), respectively.
#' @export
#' @references 
#' Pfaffl MW, Tichopad A, Prgomet C, Neuvians TP. Biotechnol Lett (2004) <doi: 10.1023/B:BILE.0000019559.84305.47>
#' @examples
#' FIBct
#' bestKeeper(FIBct)
bestKeeper=function(expression,ctVal=TRUE)  
{
  if (!ctVal){expression=log2(expression)}
  rz=list(CP.statistics=cpSta(expression,ctVal),pair.Wise.cor=pearsonCor(expression,ctVal),HKG.vs.BestKeeper=bki(expression,ctVal))
  return(rz)
}