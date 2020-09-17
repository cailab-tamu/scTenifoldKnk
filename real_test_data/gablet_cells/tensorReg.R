library(rTensor)
set.seed(1)
X1 <- matrix(runif(10*10), nrow = 10)
X1 <- (X1 + t(X1))/2
image(X1)

# X2 <- X1 * 1.1
# X3 <- X1 * 1.2
# X4 <- X1 * 1.3
# X5 <- X1 * 1.4
# X6 <- X1 * 1.5
# X7 <- X1 * 1.6
# X8 <- X1 * 1.7
# X9 <- X1 * 1.8
# X0 <- X1 * 1.9

X3 <- X5 <- X0 <- X1
X3[1:5,1:5] <- X1[1:5,1:5] * 2
X5[1:5,1:5] <- X1[1:5,1:5] * 2.5
X0[1:5,1:5] <- X1[1:5,1:5] * 3

XC <- matrix(0,10,10)
XC[1:5,1:5] <- 1
image(XC)
XC <- c(XC)

Tensor <- array(dim = c(10,10,4))
Tensor[,,1] <- X1
Tensor[,,2] <- X3
Tensor[,,3] <- X5
Tensor[,,4] <- X0
Tensor

library(rTensor)
Tensor <- as.tensor(Tensor)
set.seed(1)
oTensor <- cp(Tensor, num_components = 20)
str(oTensor$U)

plot(c(1,2,2.5,3),rowMeans(oTensor$U[[3]]))
abline(lm(rowMeans(oTensor$U[[3]])~c(1,2,2.5,3)))

T <- c(1,2,2.5,3)#rowMeans(oTensor$U[[3]])
A <- (c(oTensor$est@data[,,1]))
B <- (c(oTensor$est@data[,,2]))
C <- (c(oTensor$est@data[,,3]))
D <- (c(oTensor$est@data[,,4]))

DF <- data.frame(W=c(A,B,C,D), P = as.factor(rep(c(1:100), 4)), T= (rep(T, each = 100)))
COEF <- (coefficients(lm(W~0+poly(T,2)+P,DF)))
COEF <- COEF[-1:2]
image(abs(matrix(COEF,nrow = 10, byrow = TRUE)))
ABSCOEF <- abs(matrix(COEF,nrow = 10, byrow = TRUE))
ABSCOEF <- ABSCOEF/max(ABSCOEF)
ComplexHeatmap::Heatmap(ABSCOEF, row_order = 1:10, column_order = 1:10)
ComplexHeatmap::Heatmap(ABSCOEF)

boxplot(COEF ~ as.factor(XC))

