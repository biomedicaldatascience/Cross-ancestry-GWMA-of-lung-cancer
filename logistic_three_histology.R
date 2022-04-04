
#the logistic regression model
#input dosage was from imputed data in 0-2 format, ~570 chunks
#covariates were included: PC1-PC10, STUDY as factor in model
#CEU, CHB, YRI were stratified
#histology was stratified: adeno, squam and sclc

args = commandArgs(trailingOnly=TRUE)
chr = as.character(args[1])
chunk=as.character(args[2])

pheno=read.table("/bmds/data/Users/xxiao/list-files/70639_Disease_20PCs_10_23_2019",sep=',',head=F, na.strings=c("NA"))
names(pheno)=c("SAMPLE","D",paste("PC",1:20,sep=""),"STUDY","CEU","CHB","YRI")
pheno$STUDY=factor(pheno$STUDY)
pheno$INDEX=1:dim(pheno)[1]
hist=read.table("/bmds/data/Users/xxiao/list-files/9.studies.histology",head=F)
names(hist)=c("SAMPLE","HIST")
pheno_hist=merge(pheno,hist)
pheno_hist=pheno_hist[order(pheno_hist$INDEX),]

hist = c("adeno","squam","sclc")
temp1 = paste("/bmds/data/Users/xxiao/merge/plink/merge.",chr,".",chunk,".plink",sep="")

for (target in hist){
temp2 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy-",target,"/CEU.",chr,".",chunk,sep="")
temp3 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy-",target,"/CHB.",chr,".",chunk,sep="")
temp4 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy-",target,"/YRI.",chr,".",chunk,sep="")
if(file.exists(temp2)) {file.remove(temp2)}
if(file.exists(temp3)) {file.remove(temp3)}
if(file.exists(temp4)) {file.remove(temp4)}
}

temp51 = pheno_hist[!is.na(pheno_hist$CEU),]
temp52 = pheno_hist[!is.na(pheno_hist$CHB),]
temp53 = pheno_hist[!is.na(pheno_hist$YRI),]

xxiao = file(temp1, "rt")
while (TRUE) {
rl = readLines(xxiao, n=1)
if(length(rl)==0) {
    close(xxiao)
    break
    }

else{

    marker = strsplit(rl," ")[[1]][1]
    dose= as.numeric(c(strsplit(rl," ")[[1]][-1:-3]))

for (target in hist){

temp2 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy-",target,"/CEU.",chr,".",chunk,sep="")
temp3 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy-",target,"/CHB.",chr,".",chunk,sep="")
temp4 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy-",target,"/YRI.",chr,".",chunk,sep="")

CEU  =temp51[temp51$HIST==target | temp51$D==0,]
CEU_1=temp51[temp51$HIST==target,]
CEU_0=temp51[temp51$D==0,]

CHB  =temp52[temp52$HIST==target | temp52$D==0,]
CHB_1=temp52[temp52$HIST==target,]
CHB_0=temp52[temp52$D==0,]

YRI  =temp53[temp53$HIST==target | temp53$D==0,]
YRI_1=temp53[temp53$HIST==target,]
YRI_0=temp53[temp53$D==0,]

dose_CEU   = dose[CEU$INDEX]
dose_CEU_1 = dose[CEU_1$INDEX]
dose_CEU_0 = dose[CEU_0$INDEX]

dose_CHB   = dose[CHB$INDEX]
dose_CHB_1 = dose[CHB_1$INDEX]
dose_CHB_0 = dose[CHB_0$INDEX]

dose_YRI   = dose[YRI$INDEX]
dose_YRI_1 = dose[YRI_1$INDEX]
dose_YRI_0 = dose[YRI_0$INDEX]

AF_CEU  = sum(dose_CEU)/(2*length(dose_CEU))
AF_CEU_1= sum(dose_CEU_1)/(2*length(dose_CEU_1))
AF_CEU_0= sum(dose_CEU_0)/(2*length(dose_CEU_0))

AF_CHB = sum(dose_CHB)/(2*length(dose_CHB))
AF_CHB_1 = sum(dose_CHB_1)/(2*length(dose_CHB_1))
AF_CHB_0 = sum(dose_CHB_0)/(2*length(dose_CHB_0))

AF_YRI   = sum(dose_YRI)/(2*length(dose_YRI))
AF_YRI_1 = sum(dose_YRI_1)/(2*length(dose_YRI_1))
AF_YRI_0 = sum(dose_YRI_0)/(2*length(dose_YRI_0))

    if(var(dose_CEU) !=0){
        model1 = glm(CEU$D ~ dose_CEU+CEU$PC1+CEU$PC2+CEU$PC3+CEU$PC4+CEU$PC5+CEU$PC6+CEU$PC7+CEU$PC8+CEU$PC9+CEU$PC10+CEU$STUDY, family="binomial")
        
        result = unlist(summary(model1)[["coefficients"]][2,])
        result[5] = exp(result[1])
        cat(c(marker,AF_CEU,AF_CEU_1,AF_CEU_0,result[c(5,2,3,4)]),file=temp2,append=T,"\n")
        }
   else cat(c(marker,AF_CEU,AF_CEU_1,AF_CEU_0,NA,NA,NA,NA),file=temp2,append=T,"\n")

   
    if(var(dose_CHB) !=0){
        model1 = glm(CHB$D ~ dose_CHB+CHB$PC1+CHB$PC2+CHB$PC3+CHB$PC4+CHB$PC5+CHB$PC6+CHB$PC7+CHB$PC8+CHB$PC9+CHB$PC10+CHB$STUDY, family="binomial")
        
        result = unlist(summary(model1)[["coefficients"]][2,])
        result[5] = exp(result[1])
        cat(c(marker,AF_CHB,AF_CHB_1,AF_CHB_0,result[c(5,2,3,4)]),file=temp3,append=T,"\n")
        }
   else cat(c(marker,AF_CHB,AF_CHB_1,AF_CHB_0,NA,NA,NA,NA),file=temp3,append=T,"\n")


    if(var(dose_YRI) !=0){
        model1 = glm(YRI$D ~ dose_YRI+YRI$PC1+YRI$PC2+YRI$PC3+YRI$PC4+YRI$PC5+YRI$PC6+YRI$PC7+YRI$PC8+YRI$PC9+YRI$PC10+YRI$STUDY, family="binomial")
        
        result = unlist(summary(model1)[["coefficients"]][2,])
        result[5] = exp(result[1])
        cat(c(marker,AF_YRI,AF_YRI_1,AF_YRI_0,result[c(5,2,3,4)]),file=temp4,append=T,"\n")
        }
   else cat(c(marker,AF_YRI,AF_YRI_1,AF_YRI_0,NA,NA,NA,NA),file=temp4,append=T,"\n")


    }
  }
}

