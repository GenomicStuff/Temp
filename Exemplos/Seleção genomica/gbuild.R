#require(reshape); require(lqmm);library(MASS)

abuild<-function(ped){
  ped[,1]<-as.factor(ped[,1])
  ped[,2]<-as.factor(ped[,2])
  ped[,3]<-as.factor(ped[,3])
  
  parent<-unique(c(levels(ped[,2]),levels(ped[,3])))
  id<-ped[,1]
  
  A<-diag(length(id))
  rownames(A)<-colnames(A)<-id
  for(i in 1:length(parent)){
    id.p1<-droplevels(ped[ped[,2]%in%parent[i],1])
    id.p2<-droplevels(ped[ped[,3]%in%parent[i],1])
    halfsib<-unique(c(levels(id.p1),levels(id.p2)))
    
    A[id%in%halfsib,id%in%halfsib]<-A[id%in%halfsib,id%in%halfsib]+0.25
  }; 
  diag(A)<-1
  return(A)
}

gbuild <- function(M){
M <- as.matrix(M)
p <- NULL
for(i in 1:ncol(M)){
p <- c(p,sum(M[,i])/(2*length(M[,i])))
}

P <- matrix(p,nrow(M),ncol(M),byrow=TRUE)
P <- 2*P; W <- M-P #Parametrização (marca-2p)
#W<-M-1

G <- matrix(W%*%t(W),nrow(M),nrow(M))
G <- G/sum(2*p*(1-p))
rownames(G) <- rownames(M)
colnames(G) <- rownames(M)

#G <- G+(1-mean(diag(G)))

return(G)
}


dbuild <- function(M){
M <- as.matrix(M)
p <- NULL
for(i in 1:ncol(M)){
p <- c(p,sum(M[,i])/(2*length(M[,i])))
}

M<-rescale(abs(M-1),c(1,0))

P <- matrix(p*(1-p),nrow(M),ncol(M),byrow=TRUE)
P <- 2*P; H <- M-P #Parametrização (marca-2pq)

D <- matrix(H%*%t(H),nrow(M),nrow(M))
#D <- D/sum((2*p*(1-p)^2)) #errado!
D <- D/(2*sum(p*(1-p)*(1-2*p*(1-p))))

rownames(D) <- rownames(M)
colnames(D) <- rownames(M)

#D <- D+(1-mean(diag(D)))

return(D)
}



pbuild <- function(M){
M <- as.matrix(M)
p <- NULL
for(i in 1:ncol(M)){
p <- c(p,sum(M[,i])/(2*length(M[,i])))
}
return(p)
}

library("scales")
integer_breaks <- function(n = 5, ...) {
  breaker <- pretty_breaks(n, ...)
  function(x) {
     breaks <- breaker(x)
     breaks[breaks == floor(breaks)]
  }
}
