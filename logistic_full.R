#the logistic regression model
#input dosage was from imputed data in 0-2 format, ~570 chunks
#covariates were included: PC1-PC20, STUDY as factor in model
#CEU, CHB, YRI were stratified

args = commandArgs(trailingOnly=TRUE)
chr = as.character(args[1])
chunk=as.character(args[2])

pheno=read.table("/bmds/data/Users/xxiao/list-files/70639_Disease_20PCs_10_23_2019",sep=',',head=F, na.strings=c("NA"))
names(pheno)=c("SAMPLE","D",paste("PC",1:20,sep=""),"STUDY","CEU","CHB","YRI")
pheno$STUDY=factor(pheno$STUDY)
pheno$INDEX=1:dim(pheno)[1]
CEU = pheno[!is.na(pheno$CEU),]
CHB = pheno[!is.na(pheno$CHB),]
YRI = pheno[!is.na(pheno$YRI),]

temp1 = paste("/bmds/data/Users/xxiao/merge/plink/merge.",chr,".",chunk,".plink",sep="")
temp2 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy/CEU.",chr,".",chunk,sep="")
temp3 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy/CHB.",chr,".",chunk,sep="")
temp4 = paste("/bmds/data/Users/xxiao/plink-analysis/merge-190000-dummy/YRI.",chr,".",chunk,sep="")

r1=paste("rm ",temp2,temp22,temp3,temp33,temp4,temp44)
system(r1)

xxiao = file(temp1, "rt")
while (TRUE) {
rl = readLines(xxiao, n=1)
if(length(rl)==0) {
    close(xxiao)
    break
    }

else{

    marker = strsplit(rl," ")[[1]][1]
    al1 = strsplit(rl," ")[[1]][2]
    al2 = strsplit(rl," ")[[1]][3]
    dose= as.numeric(c(strsplit(rl," ")[[1]][-1:-3]))
    dose_CEU = dose[CEU$INDEX]
    dose_CHB = dose[CHB$INDEX]
    dose_YRI = dose[YRI$INDEX]

    if(var(dose_CEU) !=0){
        model1 = glm(CEU$D ~ dose_CEU+CEU$PC1+CEU$PC2+CEU$PC3+CEU$PC4+CEU$PC5+CEU$PC6+CEU$PC7+CEU$PC8+CEU$PC9+CEU$PC10+CEU$PC11+CEU$PC12+CEU$PC13+CEU$PC14+CEU$PC15+CEU$PC16
                           +CEU$PC17+CEU$PC18+CEU$PC19+CEU$PC20+CEU$STUDY, family="binomial")
        
        result = unlist(summary(model1)[["coefficients"]][2,])
        result[5] = exp(result[1])
        cat(c(marker,al1,al2,result[c(5,2,3,4)]),file=temp2,append=T,"\n")
        }
   else cat(c(marker,al1,al2,NA,NA,NA,NA),file=temp2,append=T,"\n")

   
    if(var(dose_CHB) !=0){
        model1 = glm(CHB$D ~ dose_CHB+CHB$PC1+CHB$PC2+CHB$PC3+CHB$PC4+CHB$PC5+CHB$PC6+CHB$PC7+CHB$PC8+CHB$PC9+CHB$PC10+CHB$PC11+CHB$PC12+CHB$PC13+CHB$PC14+CHB$PC15+CHB$PC16
                           +CHB$PC17+CHB$PC18+CHB$PC19+CHB$PC20+CHB$STUDY, family="binomial")
        
        result = unlist(summary(model1)[["coefficients"]][2,])
        result[5] = exp(result[1])
        cat(c(marker,al1,al2,result[c(5,2,3,4)]),file=temp3,append=T,"\n")
        }
   else cat(c(marker,al1,al2,NA,NA,NA,NA),file=temp3,append=T,"\n")


    if(var(dose_YRI) !=0){
        model1 = glm(YRI$D ~ dose_YRI+YRI$PC1+YRI$PC2+YRI$PC3+YRI$PC4+YRI$PC5+YRI$PC6+YRI$PC7+YRI$PC8+YRI$PC9+YRI$PC10+YRI$PC11+YRI$PC12+YRI$PC13+YRI$PC14+YRI$PC15+YRI$PC16
                           +YRI$PC17+YRI$PC18+YRI$PC19+YRI$PC20+YRI$STUDY, family="binomial")
        
        result = unlist(summary(model1)[["coefficients"]][2,])
        result[5] = exp(result[1])
        cat(c(marker,al1,al2,result[c(5,2,3,4)]),file=temp4,append=T,"\n")
        }
   else cat(c(marker,al1,al2,NA,NA,NA,NA),file=temp4,append=T,"\n")

   }
}
