library("dtw")
data("aami3a")

ref <- window(aami3a, start = 0, end = 2)
test <- window(aami3a, start = 2.7, end = 5)
alignment <- dtw(test, ref)
alignment$distance
str(aami3a)
str(alignment)
plot(dtw(test, ref, k = TRUE), type = "two", off = 1, match.lty = 2, match.indices = 20)
  
alignment <- dtw(test, ref, step.pattern = asymmetric)
alignment$distance
plot(alignment$index1,alignment$index2 )
symmetric2

query <- cbind(1:10, 1)
ref <- cbind(11:15, 2)
dtw(query, ref, dist.method = "Manhattan")$distance

lm <- matrix(nrow = 6, ncol = 6, byrow = TRUE, c(
  1, 1, 2, 2, 3, 3,
  1, 1, 1, 2, 2, 2,
  3, 1, 2, 2, 3, 3,
  3, 1, 2, 1, 1, 2,
  3, 2, 1, 2, 1, 2,
  3, 3, 3, 2, 1, 2
  ))

alignment <- dtw(lm, step = asymmetric, keep = TRUE)
alignment$costMatrix
alignment$normalizedDistance
alignmentOE <- dtw(lm, step = asymmetric, keep = TRUE, open.end = TRUE)
alignmentOE$normalizedDistance

lcm <- alignment$localCostMatrix
image(x = 1:nrow(lcm), y = 1:ncol(lcm), lcm)
text(row(lcm), col(lcm), label = lcm)
lines(alignment$index1, alignment$index2)
ccm <- alignment$costMatrix
image(x = 1:nrow(ccm), y = 1:ncol(ccm), ccm)
text(row(ccm), col(ccm), label = ccm)
lines(alignment$index1, alignment$index2)

#################################
## A noisy sine wave as query
idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;


## A cosine is for reference; sin and cos are offset by 25 samples
reference<-cos(idx)
plot(reference); lines(query,col="blue");


## Find the best match
alignment<-dtw(query,reference,keep = TRUE);

## Display the mapping, AKA warping function - may be multiple-valued
## Equivalent to: plot(alignment,type="alignment")
plot(alignment$index1,alignment$index2,main="Warping function");

## Confirm: 25 samples off-diagonal alignment
lines(1:100-25,col="red")
alignment$costMatrix

####################################################


x <- c(3,2,4,0)
y <- c(4,1,3,1)
h0 <- -1
h1 <- 2

theta <- c(0, 1)
hofx <- function(x) theta[1] + theta[2] * x
Joftheta <- function(theta, x, y) sum((hofx(x) - y)^2)/(2*length(x))

Joftheta(theta, x, y)

hoftheta(6)

x1 - x2 = -2
2*x1  -x2 = 3

1, -1
2, -1

m <- matrix(1:4, 2)
m
prop.table(m, 1)


x1 <- c( 89,72,94, 69)
x2 <- x1^2

mean_normal <-function(x) (x- mean(x))/(max(x)-min(x))
x1 <- mean_normal(x1)
x2 <- mean_normal(x2)
x <- cbind(1, x1, x2)
y <- c(96, 74, 87, 78)
alpha= 0.3
m = length(y)
theta <- c(1,1,1)
(theta <- theta - alpha/m * t(x) %*% (x%*% theta - y))
(cost_fun <- 1/(2*m) * t(x%*% theta- y) %*% (x%*%theta- y))




sin(2*pi*a*t + theta) = sin(2*pi*a*t+ theta)...

2*pi*a*t = 2*pi*n, a*t = n.... n=1 --> T = 1/a,, 주기, a : frequency 1/T

time period 2*pi...