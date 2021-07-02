library(rrBLUP); library(lattice)

gen <- read.table("genotype.txt",header=TRUE)
phe <- read.table("phenotype.txt",header=TRUE)
map <- read.table("map.txt",header=TRUE)
map <- within(map, marker <- factor(marker, levels= map$marker))

head(phe)

## DATA INFORMATION
## 
## Machos: 5; Femeas:20
## Fenotipagem: 1000; Genotipagem: 2025 indiv?duos
## Ne = 16
## h2a = 0.47 p/ o fenótipo p530
## Acuracia de sele??o fenot?pica = 0.73
## Comprimento do genoma = 500cM
## Acur?cia esperada na GS = 0.67

M <- gen[as.character(phe$ind),]; #correspondecia entre fenotipagem e genotipagem
M[1:15,1:15] #mostrando uma parte da matriz M (também chamada de matriz Z)

poptr <- rownames(M)[201:1000] #a população de treinamento
popvl <- rownames(M)[1:200] #a população de validão

g.tr <- phe[phe$ind%in%poptr,]$p530 #valores da população de treinamento
g.vl <- phe[phe$ind%in%popvl,]$p530 #valores da população de validação

train <- mixed.solve(g.tr,M[poptr,]) #rodando RRBLUP

m.cat <- train$u #cat?logo de marcas

xyplot(abs(m.cat)~map$marker, groups=map$chrom,scale="free") #Manhattan plot

g.hat <- as.matrix(M[popvl,])%*%m.cat + mean(g.tr) #valores genot?picos estimados

cor(g.hat,g.vl) #"acur?cia" realizada na GS

plot(g.hat,g.vl,xlim=c(10,45),ylim=c(10,45))
abline(0,1,col="blue")
abline(lm(g.vl~g.hat), col="red", lwd=3, lty=2)

##Eficic?ncia do m?todo utilizado (RRBLUP)
cor(g.hat,g.vl)/0.67 

##Fen?mica X Gen?mica
0.73/4
cor(g.hat,g.vl)/1
