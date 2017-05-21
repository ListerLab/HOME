args<-commandArgs(TRUE)
options(warn = -1)
file1=args
df1=read.table(args,header=TRUE)

total_reads=df1[ , grepl( "h_" , names( df1 ) ) ]
mc_read=df1[ , grepl( "mc_" , names( df1 ) ) ]
prop_table=mc_read/total_reads
workingmat = cbind(prop_table,total_reads)

label1=length(df1[ , grepl( "mc_cont" , names( df1 ) ) ])
label2=length(df1[ , grepl( "mc_case" , names( df1 ) ) ])
DX=numeric(length(mc_read))
for (i in 1:label1 ) { DX[i]=1}
for (i in (label1+1):(label1+label2) ) { DX[i]=0}
DX=as.matrix(DX)

calc_pvals <- function(combmat, DX) 
  {
  tmp_w2 = (1 / (1+exp((-combmat[(length(total_reads)+1):(length(total_reads)+length(prop_table))]))))
  weights1= (tmp_w2 - .5) / .5
  
  data2=cbind(DX,combmat[1:length(prop_table)],weights1)
  data2 = data.frame(data2)
  colnames(data2) = c( 'DX', 'props','weights1')
  glm_noW = glm(DX ~ props, family=binomial(logit), data=data2,weights = weights1)
  res = summary(glm_noW)
  if(glm_noW$converged == FALSE) {
    p_value = (1-(1 - pchisq(round(res$null.deviance-res$deviance,6), df=res$df.null-res$df.residual)))}
  else{p_value = (1 - pchisq(round(res$null.deviance-res$deviance,6), df=res$df.null-res$df.residual))}
  return(p_value)
}
p_value = apply(workingmat, 1, calc_pvals, DX)
results=cbind(df1,p_value)
write.table(results,file=file1,quote = FALSE, sep = "\t",row.names = FALSE)

