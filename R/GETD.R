#' Generalized Extended Triangular Designs
#' @param n It is a natural number such that n >= 2m ; m >= 2
#' @param m It is a natural number such that m >= 2
#'@description This package contains functions named GETD() for generating a m-associate (m >= 2) class PBIB designs and its parameters based on Generalized Triangular (GT) Association Scheme. It also calculates the Information matrix, Average variance factor and canonical efficiency factor of the generated design.
#' @return This function for generates an m-associate (m >= 2) class PBIB design and its parameters, based on Generalized Triangular (GT) Association Scheme, along with its Information matrix, Average variance factor and canonical efficiency factor.
#' @importFrom utils combn
#' @examples
#'library(GETdesigns)
#'GETD(6,2)
#'@references
#'
#'1) R.C. Bose, K.R. Nair (1939)< https://www.jstor.org/stable/40383923>" Partially balanced incomplete block designs ".
#'
#'2) R.C. Bose, T. Shimamoto (1952)<doi:10.1080/01621459.1952.10501161> "Classification and analysis of partially balanced incomplete block designs with two associate classes".
#' @export
GETD=function(n,m){
if(n>=2*m && m>=2){
######Inputs needed
k=n-m+1
################finding the m-tuple matrix
matrix1<-matrix(combn(1:n,m),ncol=m,byrow=T)
############## finding (m-1)-tuple matrix
matrix2<-matrix(combn(1:n,(m-1)),ncol=m-1,byrow=T)
################ function to find the occurances of (m-1)-tuples in m-tuple matrix
pair<-function(mat,v){
  k<-c()
  a<-v
  for(i in 1:nrow(mat)){
    #i=2
    c<-c(mat[i,])
    b<-setdiff(c,a)
    #b1<-b[b>0]
    if((length(c)-length(b)-length(a))==0){
      k<-c(k,i)
    }else{
      #k<-c(k,0)
    }
  }
  return(k)
}
############## Final design generation
final<-matrix(,nrow=0,ncol=k)
for(cc in 1:nrow(matrix2)){
  final<-rbind(final,pair(matrix1,c(matrix2[cc,])))
}
message("Triangular Design")
cat("\n")
print(final)
cat("\n")
############PARAMETERS
message("Parameters Of The Triangular Design")
cat(c("\n","Number Of Treatments (v) =",max(final),"\n","Number Of Blocks (b) =",nrow(final),"\n","Number Of Replication (r) =",length(which(final==max(final))),"\n",
      "Block Size =",ncol(final)))

##############R(replication) matrix
rep_r<-length(final[which(final==final[1,1])])
R_matrix<-diag(rep_r,max(final))
#################

################(block size) matrix
block_size<-ncol(final)
K_matrix<-diag(block_size,nrow(final))
#############incidence matrix(N-v*b)
t.incidence<-function(bibd){
  bibd
  b=nrow(bibd)
  kk=ncol(bibd)
  v<-max(bibd)
  cc<-c(1:v)

  incident<-matrix(0,nrow=v, ncol=b)
  vv=1
  while(vv<=v){
    ###########raw position of a element
    x<-which(bibd %in% c(cc[vv]))

    #########################position identify only
    k=1
    while(k<=length(x)){
      if(x[k]%%b!=0){
        x[k]<-x[k]%%b
      }else{
        x[k]<-b
      }

      k=k+1
    }
    ################

    ss=1
    while(ss<=length(x)){
      incident[vv,x[ss]]<-1
      ss=ss+1
    }
    vv=vv+1
  }
  N_prime<-t(incident)
}
N_matrix<-t(t.incidence(final))
#########
for(lm in 1:m){
  if(lm==1){
    cat(c("\n","lamda",lm,"=",1))
  }else{
    cat(c("\n","lamda",lm,"=",0))
  }
}
cat("\n")
####################C matrix
C_matrix<-R_matrix-(N_matrix%*%MASS::ginv(K_matrix)%*%t(N_matrix))
cat("\n")
message("C Matrix")
cat("\n")
print(C_matrix)
cat("\n")
############
########P matrix generation
v=max(final)
p_matrix<-matrix(0,nrow=choose(v,2),ncol=max(final))
#########
elepos<-t(combn(v,2))
#####
for(i in 1:nrow(p_matrix)){
  p_matrix[i,(elepos[i,])]<-c(1,-1)
}
########## Variance covariance part
variances<-(p_matrix)%*%MASS::ginv(C_matrix)
variances<-variances%*%t(p_matrix)
#######variance part
var<-diag(variances)
########## avg variance
Avg_var<-mean(var)
message("Average Variance Factor")
cat("\n")
print(Avg_var)
cat("\n")
#####cef
CEF<-(rep_r)/2*Avg_var
#######canonical efficiency factor
eigen<-c(eigen(C_matrix)$values)

##############
p<-c()
for(i in eigen){
  if(i>(10^(-6))){
    p<-c(p,i)
  }
}
#########harmonic mean
harmonic_mean<-1/(mean(1/p))
####canonical efficiency factor (need to change replication)
cannonical_efficiency<-(1/(rep_r))*harmonic_mean
############
message("Cannonical Efficiency Factor")
cat("\n")
print(cannonical_efficiency)
  }else{
    message("Please enter the correct value, where n >= 2*m and m>=2")
  }

}
####################################


