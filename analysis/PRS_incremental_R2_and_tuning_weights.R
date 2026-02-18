args=commandArgs(TRUE)
i=as.numeric(args[1])
library(data.table)
tuned.prs1 <- as.data.frame(fread(paste0("/PRS/validate/ukbEAS.validate.PRScsx_EAS.",i,".sscore"),header=T))
tuned.prs2 <- as.data.frame(fread(paste0("/PRS/validate/ukbEAS.validate.PRScsx_EUR.",i,".sscore"),header=T))
names(tuned.prs1)[ncol(tuned.prs1)] <- "SCORE"
names(tuned.prs1)[1] <- "FID"
names(tuned.prs2)[ncol(tuned.prs2)] <- "SCORE"
names(tuned.prs2)[1] <- "FID"
tuned.prs <- merge(tuned.prs1,tuned.prs2,by=c("FID","IID"))
tuned.prs$SCORE.x <- scale(tuned.prs$SCORE.x)
tuned.prs$SCORE.y <- scale(tuned.prs$SCORE.y)
tuned.prs <- as.data.frame(tuned.prs)[,colnames(tuned.prs)%in%c("FID","IID","SCORE.x","SCORE.y")]
tuned.ID <- as.data.frame(fread(paste0("/PRS/validateID.",i,".list"),header=F,col.names=c("FID","IID")))
pheno <- as.data.frame(fread("/PRS/ukbEAS_pheno.txt",header=T))
pheno$EA <- scale(pheno$EA)
tuned.Pheno <- merge(tuned.ID,pheno,by=c("FID","IID"))
tuned.Info <- merge(tuned.prs,tuned.Pheno,by=c("FID","IID"))
model1 <- lm(tuned.Info$EA ~., data=tuned.Info[,!colnames(tuned.Info) %in% c("FID","IID","EA")])
coeffX <- summary(model1)$coefficients[2,1]
coeffY <- summary(model1)$coefficients[3,1]

test.prs1 <- as.data.frame(fread(paste0("/PRS/test/ukbEAS.test.PRScsx_EAS.",i,".sscore"),header=T))
test.prs2 <- as.data.frame(fread(paste0("/PRS/test/ukbEAS.test.PRScsx_EUR.",i,".sscore"),header=T))
names(test.prs1)[ncol(test.prs1)] <- "SCORE"
names(test.prs1)[1] <- "FID"
names(test.prs2)[ncol(test.prs2)] <- "SCORE"
names(test.prs2)[1] <- "FID"
test.prs <- merge(test.prs1,test.prs2,by=c("FID","IID"))
test.prs$SCORE.x <- scale(test.prs$SCORE.x)
test.prs$SCORE.y <- scale(test.prs$SCORE.y)
test.prs <- as.data.frame(test.prs)[,colnames(test.prs)%in%c("FID","IID","SCORE.x","SCORE.y")]
test.prs$PRS <- coeffX*test.prs$SCORE.x + coeffY*test.prs$SCORE.y
test.ID <- as.data.frame(fread(paste0("/PRS/testID.",i,".list"),header=F,col.names=c("FID","IID")))
test.Pheno <- merge(test.ID,pheno,by=c("FID","IID"))
test.Info <- merge(test.prs,test.Pheno,by=c("FID","IID"))

model_null <- lm(test.Info$EA ~., data=test.Info[,!colnames(test.Info)%in%c("FID","IID","SCORE.x","SCORE.y","EA","PRS")])
model_test <- lm(test.Info$EA ~., data=test.Info[,!colnames(test.Info)%in%c("FID","IID","SCORE.x","SCORE.y","EA")])
model_prs1 <- lm(test.Info$EA ~., data=test.Info[,!colnames(test.Info)%in%c("FID","IID","SCORE.y","EA","PRS")])
model_prs2 <- lm(test.Info$EA ~., data=test.Info[,!colnames(test.Info)%in%c("FID","IID","SCORE.x","EA","PRS")])
r2_null <- summary(model_null)$r.squared
r2_test <- summary(model_test)$r.squared
r2_prs1 <- summary(model_prs1)$r.squared
r2_prs2 <- summary(model_prs2)$r.squared
result_test <- data.frame(R2=r2_test-r2_null,SE=NA,weight_EAS=coeffX,weight_EUR=coeffY,iter=i)
library(boot)
lm_model <- function(data,indices) {
  null <- lm(data$EA ~., data = data[indices,!colnames(data)%in%c("FID","IID","SCORE.x","SCORE.y","PRS","EA")])
  model <- lm(data$EA ~., data=data[indices,!colnames(data)%in%c("FID","IID","SCORE.x","SCORE.y","EA")])
  return(summary(model)$r.squared-summary(null)$r.squared)
}
bootstrap_results <- boot(test.Info, lm_model, R = 1000)
result_test$SE <- sd(bootstrap_results$t)
fwrite(result_test,paste0("/PRS/ukbEAS.test.PRScsx.",i,".R2.result"),quote=F,sep="\t",col.names=T,row.names=F,na="NA")

