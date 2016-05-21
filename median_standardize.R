standardize <- function(z) {
  rowmed <- apply(z, 1, median)
  rowmad <- apply(z, 1, mad)  # median absolute deviation
  rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
  rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
  return(rv)
}

>data<-data.frame(data)
>colnames(data)<-c("sample1","sample2","sample3","sample4")
>rownames(data)<-c("gene1","gene2","gene3")
>standardize(data)
sample1   sample2    sample3   sample4
gene1 -1.011736 -0.3372454 0.3372454 1.011736
gene2 -1.011736 -0.3372454 0.3372454 1.011736
gene3 -1.011736 -0.3372454 0.3372454 1.01173