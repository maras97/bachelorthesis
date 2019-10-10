library(ROCR)
library(glmnet)

get.import.coef <- function(glmnet.fit, lambda){
  glmcoef<-coef(glmnet.fit,lambda)
  coef.increase<-dimnames(glmcoef[glmcoef[,1]>0,0])[[1]]
  coef.decrease<-dimnames(glmcoef[glmcoef[,1]<0,0])[[1]]
  #get ordered list of variables as they appear at smallest lambda
  allnames<-names(coef(glmnet.fit)[,
                                   ncol(coef(glmnet.fit))][order(coef(glmnet.fit)[,
                                                                                  ncol(coef(glmnet.fit))],decreasing=TRUE)])
  #remove intercept
  allnames<-setdiff(allnames,allnames[grep("Intercept",allnames)])
  
  #assign colors
  cols<-rep("gray",length(allnames))
  cols[allnames %in% coef.increase]<- "deeppink4" #"darkseagreen3"     
  cols[allnames %in% coef.decrease]<- "dodgerblue"
  return(cols)
}



## DATA PREPARATION
##
# reading motifcounter results for high- and low- scored sequences from mESC with detailed information on scores for every sequence and motif cluster
matX <- readRDS("../RESULTS_motifcounter/motifcounter_StrongWeakESC_detailed_clusteredMotifs.rds") # 514 rows, 79 cols
Y <- c(rep(1,257), rep(0,257)) #1: high-scored (active enhancer), 0:low-scored (inactive enhancer)
set.seed(25) # set seed for random number generator to get reproducable results

# dataframe w/ predictor variables & a response variable
# 1st col : response variable; # cols 2-80 = predictor variables
df <- (data.matrix(cbind(Y,matX))) 

# Create training subset for model development & testing set for model performance testing
inTrain <- createDataPartition(df[,1], p = .75, list = FALSE)
Train <- df[ inTrain, ]
Test <- df[ -inTrain, ]
#Train <- readRDS("../classifier/TrainSet.rds")
#Test <- readRDS("../classifier/TestSet.rds")


## PARAMETER SELECTION 
##
## NOTE: cross-validation picks k-fols randomly (independent of seed) so that it produces slighty different results each time: to reproduce results save obtained objects in RDS-files  

## plot a roc curve for different alpha-parameter values 
alpha <- c(0,0.5,0.9,1)
colour <- c("red","orange", "green","blue")
auc_all <- c()
pdf("../classifier/ROC_varAlpha.pdf")
for (a in 1:length(alpha)){
  # get an optimal lambda for the given alpha value by cross-validation and then test the model on the testing set 
  lasso.model <- cv.glmnet(x = Train[,2:80], y = Train[,1], 
                           family = 'binomial', type.measure = 'auc', alpha = alpha[a], intercept = F)
  lasso.prob <- predict(lasso.model,type="response", 
                        newx = Test[,2:80], s = lasso.model$lambda.1se)
  pred <- prediction(lasso.prob, Test[,1])
  # calculate probabilities for TPR/FPR for predictions
  perf <- performance(pred,"tpr","fpr")
  # get AUC for model
  q <- performance(pred,"auc") 
  auc <- q@y.values[[1]]
  auc_all <- append(auc_all,round(auc, digits = 2))
  c <- colour[a]
  plot(perf,colorize=F, col=c ,lwd=2) # plot ROC curve
  par(new=TRUE)
}
lines(c(0,1),c(0,1),col = "gray", lty = 4)
title(paste("AUC alpha=0.0:", format(auc_all[1],nsmall=2), "\nAUC alpha=0.5:", format(auc_all[2],nsmall=2), "\nAUC alpha=0.9:", format(auc_all[3],nsmall=2),"\nAUC alpha=1.0:", format(auc_all[4],nsmall=2)), font.main = 1, cex.main = 0.8)
legend("bottomright",legend=c("alpha = 0","alpha = 0.5","alpha = 0.9","alpha = 1"), col = colour, lty = 1, cex = 0.8, lwd = 2)
dev.off()


## calculate AUC for different alpha- and lambda-values and plot in a heatmap to obtain an overview of model performance dependent on alpha & lambda 
alpha <- rev(c(seq(0,0.5,0.1),0.75,1))
lambda <- seq(0,0.4,0.01)
max <- c()
auc_all <- matrix(nrow=length(alpha),ncol=length(lambda)) # rows: alpha, cols: lambda 
for (a in 1:length(alpha)){
  for (l in 1:length(lambda)){
    # train model with the given alpha- and lambda-value and test model with the testing set
    lasso.model <- cv.glmnet(x = Train[,2:80], y = Train[,1], 
                             family = 'binomial', type.measure = 'auc', alpha = alpha[a], intercept = F)
    lasso.prob <- predict(lasso.model,type="response", 
                          newx = Test[,2:80], s = lambda[l])
    pred <- prediction(lasso.prob, Test[,1])
    # calculate probabilities for TPR/FPR for predictions
    perf <- performance(pred,"tpr","fpr")
    # get AUC for model
    q <- performance(pred,"auc") 
    auc <- q@y.values[[1]]
    print(auc_all)
    auc_all[a,l] <- auc
    if (auc==max(auc_all, na.rm=T)){
      max[1] <- alpha[a]
      max[2] <- lambda[l]
    }
  }
}
rownames(auc_all) <- format(alpha)
colnames(auc_all) <- format(lambda)
# plot obtained matrix with AUCs in a heatmap with colums = lambda and rows = alpha
pal <- c("oldlace","lightgoldenrodyellow", "cyan3" ,"darkblue")
col_fun = colorRamp2(c(min(auc_all),(min(auc_all)+ max(auc_all))/2,(min(auc_all)+ max(auc_all)*2)/3, max(auc_all)),col = pal)
pdf("../classifier/HeatmapParamters_varAlphaLambda.pdf")
h <- Heatmap(auc_all, cluster_rows = F, cluster_columns = F, row_title_side = "left", row_names_side = "left",
             col = col_fun, name = "AUC", column_title = "lambda : Amount of Regularization", column_names_gp = gpar(fontsize = 7), 
             column_title_side = "bottom",
             row_title = "alpha : Mixing Percentage \n Ridge <-----> Lasso")
draw(h, column_title =  paste("maximal AUC = ", round(max(auc_all), digits= 2), " for alpha = ", max[1], " & lambda = ", max[2], sep=""), 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


## TRAINING OF THE LOGISTC REGRESSION MODEL
##
## plot to show performance dependent on lambda
cvfit = cv.glmnet(x=Train[,2:80], y=Train[,1], family = "binomial", type.measure = "class", alpha=0.1, intercept = F)
#cvfit <- readRDS("../classifier/cvfit.rds")
pdf("../classifier/glmnetFit_alpha0.1.pdf")
plot(cvfit)
title(paste("alpha =",0.1), line = 3, font.main=1, cex.main = 0.95) 
dev.off()

# plot to see the development of coefficients dependent on lambda 
fit = glmnet(x=Train[,2:80], y=Train[,1], family = "binomial", alpha=0.1, intercept = F)
#saveRDS(fit, "../fit.rds")
#fit <- readRDS("../classifier/fit.rds")
colors <- get.import.coef(fit,cvfit$lambda.1se) # |coeffs|>0 for lambda.1se will be highlighted
pdf("../classifier/glmnetLambdaPlot_alpha0.1_lambda1se.pdf")
plot_glmnet(fit, label = T, s = cvfit$lambda.1se, col = colors, main = "alpha = 0.1")
dev.off()

## TESTING PERFORMANCE OF THE LOGISTC REGRESSION MODEL
##
# Apply model to testing dataset
lasso.prob <- predict(cvfit,type="response", newx = Test[,2:80], s = cvfit$lambda.1se)
pred <- prediction(lasso.prob, Test[,1])
# calculate probabilities for TPR/FPR for predictions
perf <- performance(pred,"tpr","fpr")
# get AUC for model
q <- performance(pred,"auc") # shows calculated AUC for model
auc <- q@y.values[[1]]
#  0.8168945 für lambda_1se, alpha = 0.1

# sample the random seed 100-times to split the data differently into training- and testing set to validate AUC of regression model
rand.seeds <- sample(100)
aucs <- c()
for (k in rand.seeds){
  inTrain <- createDataPartition(df[,1], p = .75, list = FALSE)
  Train <- df[ inTrain, ]
  Test <- df[ -inTrain, ]
  cvfit = cv.glmnet(x=Train[,2:80], y=Train[,1], family = "binomial", type.measure = "class", alpha=0.1, intercept = F )
  lasso.prob <- predict(cvfit,type="response", newx = Test[,2:80])
  pred <- prediction(lasso.prob, Test[,1])
  # calculate probabilities for TPR/FPR for predictions
  perf <- performance(pred,"tpr","fpr")
  # get AUC for model
  q <- performance(pred,"auc") 
  auc <- q@y.values[[1]]
  aucs <- append(aucs, auc)
}
# calcualte mean AUC from 100 times seed-sampling
mean(aucs)



## EXPORT BETA COEFFICIENTS OF THE LOGISTC REGRESSION MODEL
##
# extract beta coefficients from the fitted regression model and sort in descending order 
coefs <- coef(fit, s= cvfit$lambda.1se)
coefs_ordered <- coefs[order(abs(coefs), decreasing = T),]
coefs_ordered <- coefs_ordered[coefs_ordered!=0]
# load representative TFs for cluster and match with coefs to export in a table with rounded coefficients associated with the motif clusters
cl_id <- read.table("../motifDB/Clusters_1TF.tab", col.names = c("motif", "matrix", "TF"))
tf_id <- match(names(coefs_ordered),cl_id[,1])
names_cluster <- cl_id[tf_id,3]
clShort <- paste(names(coefs_ordered),names_cluster,sep=":")
matCoefs <- cbind(clShort, round(coefs_ordered, digits=5))
write.csv(matCoefs, "../classifier/glmnetLambda1se_coefficientsRounded.csv" ,quote = F,  row.names = F)





