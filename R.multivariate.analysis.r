### linear regressions ###
wine <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data",
                   sep=",")

head(wine)


plot(wine$V7, wine$V8)

### linear regressions ###
### Y = a * X + b ###
### Y --> target variable
### X --> predictor variable


fit = lm(V8 ~ V7, data=wine)
summary(fit)



### V8 --> target variable
### V7 --> predictor variable


### V7, estimate ==> a coefficient
### (intercept), estimate ==> b coefficient
a = 1.37984
b = -1.13763

V8.predicted = a * wine$V7 + b

plot(V8.predicted, wine$V8) ### predicted Y vs actual Y 
abline(a=0,b=1, col="red") ### drawing diagnoal line


### predict function ###

V8.predicted = predict(fit)
plot(V8.predicted, wine$V8)
abline(a=0,b=1, col="red") #diagnoal line


### drawing regression lines over Y ~ X ###

plot (wine$V7, wine$V8)

plot(wine$V7, wine$V8)
abline(lm(V8 ~ V7, data=wine), col="red")
summary(lm(V8~V7, data=wine))
### adjusted R-square: 0.746



### multiple regressions ###
### Y = a1 * X1 + a2 * X2 + a3 * X3 + b
fit.multiple = lm(V8 ~ V2 + V3 + V4 + V5 + V6 + V7, data=wine)
AIC(fit.multiple) # 256.404024437778
# Akaike information criterion = AIC 
# Lower AIC is better in that model is simpler and better fit for the data

plot(predict(fit.multiple, wine), wine[,"V8"])
abline(a=0,b=1, col="red")


summary(fit.multiple)
### adjusted R-square: 0.7642



fit.multiple2 = lm(V8 ~ V3 + V4 + V5 + V7, data=wine)
AIC(fit.multiple2) # 253.01534015914


plot(predict(fit.multiple2, wine), wine[,"V8"])
abline(a=0,b=1, col="red")


summary(fit.multiple2)
### adjusted R-square: 0.7662



fit.multiple3 = lm(V8 ~ V3 + V5 + V7, data=wine)
summary(fit.multiple3)
AIC(fit.multiple3) # 254.57126639392


plot(predict(fit.multiple3, wine), wine[,"V8"])
abline(a=0,b=1, col="red")


summary(fit.multiple3)
### adjusted R-square: 0.7628



eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}
### Lasso regression, regression that gives a penalty to complex model, choose best features (L2 regularization)  ###
if (!require(glmnet)) install.packages("glmnet")
require(glmnet)
wine.mat = as.matrix(wine)

wine.train = wine.mat[1:100,]
wine.test  = wine.mat[101:178,]

wine.train.y = wine.train[,"V8"]
wine.train.x = wine.train[,colnames(wine.train) != "V8"]

wine.test.y = wine.test[,"V8"]
wine.test.x = wine.test[,colnames(wine.test) != "V8"]


lambdas <- 10^seq(2, -3, by = -.1)

# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(wine.train.x, wine.train.y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)

# Best 
lambda_best <- lasso_reg$lambda.min 
lambda_best

lasso_model <- glmnet(wine.train.x, wine.train.y, alpha = 1, lambda = lambda_best, standardize = TRUE)
#summary(lasso_model$beta)
#dim(wine)

predictions_train <- predict(lasso_model, s = lambda_best, newx = wine.train.x)
plot(predictions_train, wine.train.y)
abline(a=0,b=1,col="red")
eval_results(wine.train.y,  predictions_train, wine.train)


predictions_test <- predict(lasso_model, s = lambda_best, newx = wine.test.x)
plot(predictions_test, wine.test.y)
abline(a=0,b=1,col="red")
eval_results(wine.test.y,  predictions_test, wine.test)


### Eigenvector ###
###

A <- matrix(c(13, -4, 2, -4, 11, -2, 2, -2, 8), 3, 3, byrow=TRUE)
A
colnames(A) = c("A1","A2","A3")

###
### Matrix can be expressed as multiple combinations of vectors
### which is eigenvectors ###
###
### A * v = Lambda * v
###
### v = eigenvector
### Lambda = eigenvalue
###
### then why we need????
### v, eigenvector, is kind of axis of space that matrix A produces
###




ev <- eigen(A)
ev

eigen.vectors = ev$vectors
colnames(eigen.vectors) = c("E1","E2","E3")


eigen.vectors



### let's check eigen vector ###

#example.mat = cbind(A, eigen.vectors)

A

# A * v
multiplied =   A %*% eigen.vectors[,1]

# lambda * v
multiplied2 = ev$values[1] * eigen.vectors[,1]


data.frame(A_v = multiplied, lambda_v = multiplied2)


# A * v
multiplied =   A %*% eigen.vectors[,3]


# lambda * v
multiplied2 = ev$values[3] * eigen.vectors[,3]


data.frame(A_v = multiplied, lambda_v = multiplied2)




library(ggplot2)
set.seed(1)
base1 <- rnorm(20)
x1 <- base1 + rnorm(20) * 0.4
y1 <- base1 + rnorm(20) * 0.4

dat <- data.frame(id = 1:20, x1, y1)

ggplot(dat, aes(x = x1, y = y1)) + geom_point(size = 5) + 
  geom_text(aes(label = id), col = "white") + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  geom_abline(slope = -1, intercept = 0, linetype = 2) + 
  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)

### rotation by multiplications with eigen vectors

e <- eigen(cov(dat[, -1]))

rotated <- data.frame(id = dat$id, as.matrix(dat[, 2:3]) %*% e$vectors)

ggplot(rotated, aes(x = X1, y = X2)) + 
  geom_point(size = 5) + 
  geom_text(aes(label = id), col = "white") +
  geom_abline(intercept = 0, slope = 0)


### eigen decomposition --> PCA ###

pca.out = prcomp(dat[,2:3])


pca.result <- data.frame(id = dat$id, X1 = pca.out$x[,1], X2 = pca.out$x[,2])

ggplot(pca.result, aes(x = X1, y = X2)) + 
  geom_point(size = 5) + 
  geom_text(aes(label = id), col = "white") + 
  geom_abline(intercept = 0, slope = 0)

### PCA analysis ###

library(datasets)
data(iris)


### principal component analysis ###
pca.out = prcomp(iris[,1:4])
plot(pca.out$x[,1:2], col=iris$Species, pch=16)



### multidimensional scaling (MDS) ###
iris.dist = dist((iris[,1:4]), method="euclidean")
mds.out = cmdscale(iris.dist, eig = T, k = 2)


plot(mds.out$points[,1], mds.out$points[,2], col=iris$Species, pch=16, main = "MDS plot")


#### non-metric multidimensional scaling (NMDS) ####


require(vegan)
set.seed(2)

nmds.out = metaMDS(comm = iris[,1:4],k=2, trymax=100)
plot(nmds.out$points, col=iris$Species, xlab="NMDS1", ylab="NMDS2", pch=16)
stressplot(nmds.out)


### tSNE analysis ###
require(Rtsne)
set.seed(1)
tsne.out = Rtsne(iris[,1:4], check_duplicates = F)
plot(tsne.out$Y[,1:2], col=iris$Species, xlab = "tSNE1", ylab = "tSNE2", pch=16 )




### UMAP analysis ###


require(umap)

umap.out = umap(iris[,1:4])

plot(umap.out$layout[,1:2], col=iris$Species, xlab = "UMAP1", ylab = "UMAP2", pch=16)



