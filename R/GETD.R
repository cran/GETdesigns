#' Generalized Extended Triangular Designs
#' @param n It is a natural number such that n >= 2m ; m >= 2
#' @param m It is a natural number such that m >= 2
#' @param trt Provide any treatment number to know its all associates. By default it is 1.
#'@description This package contains a function named GETD() for generating m-associate (m>=2) class PBIB designs along with parameters (v, b, r, k and lambda_i, i = 1, 2,…,m) and the underlying Generalized Extended Triangular (GET) Association Scheme.
#' @return This package generates an m-associate (m >= 2) class PBIB designs under GET Association Scheme. It also calculates the Information matrix, Average variance factor, canonical efficiency factor and different treatment associates of the generated designs.
#' @importFrom utils combn
#' @examples
#'library(GETdesigns)
#'GETD(6,2,1)
#'@references
#'
#'1) R.C. Bose, K.R. Nair (1939)< https://www.jstor.org/stable/40383923>" Partially balanced incomplete block designs ".
#'
#'2) R.C. Bose, T. Shimamoto (1952)<doi:10.1080/01621459.1952.10501161> "Classification and analysis of partially balanced incomplete block designs with two associate classes".
#' @export
GETD=function(n,m,trt=1){
  options("max.print"=100000)
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

# ####to find associates
# as_list<-list()
# total<-c(1:nrow(matrix1))
# for(j in 1:nrow(matrix1)){
#   for(l in 1:(m)){
#   k=1
#   for(i in total[total!=j]){
#
#       if(length(setdiff(matrix1[j,],matrix1[i,]))==length(matrix1[j,])-l){
#         k=c(k,j)
#       }
#     }
#   }
#   as_list=append(as_list,list(k))
# }

######################
for(lm in 1:m){
  if(lm==1){
    cat(c("\n","lamda",lm,"=",1))
  }else{
    cat(c("\n","lamda",lm,"=",0))
  }
}
cat("\n")
####Parameters of Association schemes
message("Parameters of Association scheme")
cat("\n")
###Associates for a treatment
trt_list<-list()
for(l in (m-1):0){
  store_trt=c()
  select_row=matrix1[trt,]
  for(o in 1:nrow(matrix1)){
    if(length(setdiff(select_row,matrix1[o,]))==ncol(matrix1)-l){
      store_trt=c(store_trt,o)
    }
  }
  trt_list<-append(trt_list,list(store_trt))
}
#########
##associates
english<-c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth", "Ninth", "Tenth",
          "Eleventh", "Twelfth", "Thirteenth", "Fourteenth", "Fifteenth", "Sixteenth", "Seventeenth", "Eighteenth", "Nineteenth", "Twentieth",
          "Twenty-first", "Twenty-second", "Twenty-third", "Twenty-fourth", "Twenty-fifth", "Twenty-sixth", "Twenty-seventh", "Twenty-eighth", "Twenty-ninth", "Thirtieth",
          "Thirty-first", "Thirty-second", "Thirty-third", "Thirty-fourth", "Thirty-fifth", "Thirty-sixth", "Thirty-seventh", "Thirty-eighth", "Thirty-ninth", "Fortieth",
          "Forty-first", "Forty-second", "Forty-third", "Forty-fourth", "Forty-fifth", "Forty-sixth", "Forty-seventh", "Forty-eighth", "Forty-ninth", "Fiftieth",
          "Fifty-first", "Fifty-second", "Fifty-third", "Fifty-fourth", "Fifty-fifth", "Fifty-sixth", "Fifty-seventh", "Fifty-eighth", "Fifty-ninth", "Sixtieth",
          "Sixty-first", "Sixty-second", "Sixty-third", "Sixty-fourth", "Sixty-fifth", "Sixty-sixth", "Sixty-seventh", "Sixty-eighth", "Sixty-ninth", "Seventieth",
          "Seventy-first", "Seventy-second", "Seventy-third", "Seventy-fourth", "Seventy-fifth", "Seventy-sixth", "Seventy-seventh", "Seventy-eighth", "Seventy-ninth", "Eightieth",
          "Eighty-first", "Eighty-second", "Eighty-third", "Eighty-fourth", "Eighty-fifth", "Eighty-sixth", "Eighty-seventh", "Eighty-eighth", "Eighty-ninth", "Ninetieth",
          "Ninety-first", "Ninety-second", "Ninety-third", "Ninety-fourth", "Ninety-fifth", "Ninety-sixth", "Ninety-seventh", "Ninety-eighth", "Ninety-ninth", "Hundredth")
message("Association Table for Treatment ",trt)
for(i in 1:length(trt_list)){
  print(paste(english[i],"associates are:"),quote=FALSE)
  print(trt_list[[i]])
  cat("\n")
}
cat("\n")
###############ni s find
##ni s
for(i in 1:m){
  print(paste0("n_",i,"=",choose(m,i)*choose(n-m,i)),quote=FALSE)
  cat("\n")
}
cat("\n")
##########P matrices
###Pi s
list_of_p<-list()
for(i in 1:m){
  pi_jk=c()
  for(j in 1:m){
    for(k in 1:m){
      sum=c()
      for(h in 0:(m-i)){
        sum=c(sum,prod(choose(m-i,h),choose(i,m-j-h), choose(i,m-k-h),choose(n-m-i,j+k+h-m)))
      }
      pi_jk=c(pi_jk,sum(sum))
    }
  }
  list_of_p<-append(list_of_p,list(matrix(pi_jk,m,m)))
}
message("P Matrices")
print(list_of_p)
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


