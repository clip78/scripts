## calculate p-value for three circle venn diagram 
## intersection
drawtest3<-function(N,m,a,b,c){
  #N=7000 # total number of genes for selection
  #m=3 # number of overlaps in all three sets
  #a=200 # number of up-regulated genes in sample A
  #b=250 # number of up-regulated genes in sample B
  #c=300 # number of up-regualted genes in sample C
  
  #create a dataframe to store frequency table
  d<-data.frame(0:min(a,b,c),rep(0,min(a,b,c)+1))
  
  for (i in 1:10000){
    A=sample(1:N,size=a,replace=FALSE)
    B=sample(1:N,size=b,replace=FALSE)
    C=sample(1:N,size=c,replace=FALSE)
    e<-table((C %in% A)&(C %in% B))["TRUE"]
    # if there is no intersection, put 0 instead of NA
    if(!complete.cases(e)){
      e<-0
    }
    # add one to the counter for corresponding occurence
    d[e+1,2]<-d[e+1,2]+1
  }
  colnames(d)<-c("Intersect","p-value")
  # convert counts to ratios
  d[,2]<-d[,2]/10000
  # calculate p-value
  p<-sum(d[(m+1):(min(a,b,c)+1),2])
  #print(d)
  return(p)
}