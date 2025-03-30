## QUESTION-1

#Part A: Methods for density estimation 
#Univariate Density Estimation Methods: histogram, frequency polygon, average shifted histogram, and kernel density estimators: It is suitable for Tooth Growth data since univariate density estimation methods do not have assumptions other than contiunous random variable.

#Kernel Density Estimation: is a non-parametric way to estimate the probability density function of a random variable. Since the ToothGrowth data is finite data and samples are univariate independent and identically distributed sample. It is suitable using this method since we may interested in the shape of the distribution if we are not very sure our data(ToothGrowth) is normally distributed.

#Parametric probability density estimation method: We assume the population follows a certain distribution and try to estimate its parameters from the available data. The response variable of Tooth Growth data is normally distributed. So we can use parametric probability density estimation methods.

#Part B: The histogram of this sample using Sturges’ Rule

set.seed(361)
# Sturges' Rule
sturges_rule <- function(x){
  n <- length(x) #sample size 
  
  #The number of classes for a Histogram
  nclass <- ceiling(1 + log2(n)) #denominator part of the formula
  
  #The difference between each consecutive pair of elements of a vector (bin width according to Sturges' Rule)
  cwidth <- diff(range(x) / nclass)  
  
  #The optimal number of breaks (bin)
  breaks <- min(x) + cwidth * 0:nclass
  
  return(list(nclass = nclass, cwidth = cwidth, breaks = breaks))
}

#Sturges' Rule in samples which are normally distributed.(We can use since our sample is normally distributed.)
x <- ToothGrowth$len
z_sturges <- seq(min(x) - sturges_rule(x)$cwidth, max(x) + sturges_rule(x)$cwidth, 0.01)

h.sturges <- hist(x, breaks = sturges_rule(x)$breaks, prob = T, main = "Sturges' Rule Histogram")
lines(z_sturges, dnorm(z_sturges, mean(ToothGrowth$len), sd(ToothGrowth$len)), col = "pink", lwd = 5)


#Part C: The histogram of the sample using Scott’s Normal Reference Rule

scotts_rule <- function(x){
  n <- length(x) #sample size
  h <- 3.5 * sd(x) * (n^(-1/3)) #bin width according to scotts rule
  nclass <- ceiling(diff(range(x)) / h) #number of classes
  breaks <- min(x) + h * 0:nclass #breaks for histogram
  return(list(nclass = nclass, h = h, breaks = breaks))
}

x <- ToothGrowth$len
h.scott <- hist(x, breaks = scotts_rule(x)$breaks, freq = F, main = "Scotts Rule Histogram",
                ylim= c(0,0.06),xlim=c(0,40))

z_scott <- seq(min(x) - scotts_rule(x)$h, max(x) + scotts_rule(x)$h, 0.01)
lines(z_scott, dnorm(z_scott, mean(ToothGrowth$len), sd(ToothGrowth$len)), col = "pink", 
      lwd = 5)


round(h.scott$breaks, 1) #breaks
h.scott$counts # counts; number of obsv in each bar
scotts_rule(x)$h #The bin width according to Scotts' rule


#Part D: The histogram of the sample using Freedman-Diaconis Rule

FD_rule <- function(x){
  n <- length(x) #sample size
  h <- 2 * IQR(x) * n^(-1/3) #bin width according to Freedman-Diaconis Rule
  nclass <- ceiling(diff(range(x)) / h) #number of classes in histogram
  breaks <- min(x) + h * 0:nclass #breaks of histogram
  return(list(nclass = nclass, h = h, breaks = breaks))
}

x <- ToothGrowth$len
h.FD <- hist(x, breaks = FD_rule(x)$breaks, freq = FALSE, main = "Freedman-Diaconis Rule Histogram", ylim= c(0,0.06), xlim=c(0,40))
z_fd <- seq(min(x) - FD_rule(x)$h, max(x) + FD_rule(x)$h, 0.01)
lines(z_fd, dnorm(z_fd, mean(ToothGrowth$len), sd(ToothGrowth$len)), col = "pink", lwd = 5)

#Part E: The histogram of the sample using Frequency Polygon Density Estimate

freq_poly <- function(x){
  
  n <- length(x) #sample size
  h <- 2.15 * sd(x) * n^(-1/5) # bin width according to Frequency Polygon Density Estimate
  br <- pretty(x, diff(range(x)) / h) #number of classes
  brplus <- c(min(br)-h, max(br)+h) #breaks
  return(list(brplus = brplus, h = h, br = br))
  
}

x <-ToothGrowth$len

h.freq <- hist(x, breaks = freq_poly(x)$br, freq = F, main = "Frequency Polygon Density Estimate Histogram",
               xlim = freq_poly(x)$brplus)

#density est at vertices of polygon
vx <- h.freq$mids
vy <- h.freq$density
delta <- diff(vx)[1] # h after pretty is applied

k <- length(vx)
vx <- vx + delta # the bins on the ends

vx <- c(vx[1] - 2 * delta, vx[1] - delta, vx)
vy <- c(0, vy, 0)
# add the polygon to the histogram
polygon(vx, vy)

#Part F: TWO density plot by using Kernel Density Estimation method with ROBUST Bandwidth and with a KERNEL which you want to try, and find density of point x0 = 11.2.

x <- ToothGrowth$len
n <- length(x)
h <- 0.9 * min(sd(x), IQR(x)) * n ^ (-1/5)
plot(density(x))

d <- density(x, bw = h)
xnew <- seq(min(x), max(x), 0.009)
fhat <- approx(d$x, d$y, xout = xnew)
plot(fhat)


n <- length(x)
plot(density(x), main="",ylim= c(0,0.05))

abline(v = 11.2)
y <- seq(-5,40, 5)

lines(y, dnorm(y, mean(ToothGrowth$len), sd(ToothGrowth$len)), lty = 2)

#Part G: Suppose x0 = 7.2. Find which bin does x0 correspond to and what is the density for x0 according to Default histogram, Sturges’, Scott’s and Freedman Rules.

x0 <- 7.2

#Default Histogram
h.default <- hist(x, freq = FALSE,  main = "default n=25")
b <- which.min(h.default$breaks <= x0) - 1
print(c(b, h.default$density[b]))

#Sturges Rule Histogram
b <- which.min(h.sturges$breaks <= x0) - 1
print(c(b, h.sturges$density[b]))

#Scotts' Rule Histogram
b <- which.min(h.scott$breaks <= x0) - 1
print(c(b, h.scott$density[b]))

#Freedman Rule Histogram
b <- which.min(h.FD$breaks <= x0) - 1
print(c(b, h.FD$density[b]))


## QUESTION-2

#Part A
hist(PlantGrowth$weight)
shapiro.test(PlantGrowth$weight) #The data is normally distributed.
x <- PlantGrowth$weight

n <- length(x)

#The optimal bin width for the naive ASH estimate of a Normal(μ,σ^2) density is:
h <- 2.576 * sd(x) * n^(-1/5)

a <- min(x) - 0.5;b <- max(x) + 0.5

for (l in 1:5){
  
  sample_size <- c(2, 5 ,10, 50, 100)
  
  delta <- h / sample_size[l] 
  breaks <- seq(a - h, b + 2*h, delta) #breaks formula for ASH
  
  hist.ash <- hist(x, breaks = breaks, plot = F)
  nk <- hist.ash$counts
  K <- abs((1-sample_size[l]):(sample_size[l]-1))
  
  fhat <- function(x){
    i <- max(which(x > breaks))
    k <- (i - sample_size[l] + 1):(i + sample_size[l] - 1)
    vk <- nk[k]
    sum((1 - K / sample_size[l]) * vk) / (n * h)
  }
  fhat(3.6)
  
  # density can be computed at any points in range of data
  z_ash <- as.matrix(seq(a, b + h, 0.1))
  f.ash <- apply(z_ash, 1, fhat) #density estimates at midpts
  
  # plot ASH density estimate over histogram
  breaks2 <- seq(a, b + h, h)
  hist(x, breaks = breaks2, freq = F, main = "The Averaged Shifted Histogram of PlantGrowth", ylim = c(0, max(f.ash)))
  lines(z_ash, f.ash, xlab = "x")

#Part B
  h.ASH <- hist(x, breaks = breaks2, freq = F, main = "The Averaged Shifted Histogram of PlantGrowth", ylim = c(0, max(f.ash)))
  
  #When x0 = 4.1
  x0 <- 4.1
  b <- which.min(h.ASH$breaks <= x0) - 1
  print(c(b, h.ASH$density[b]))
  
  #When x1 = 3.5
  x1 <- 3.5
  b <- which.min(h.ASH$breaks <= x1) - 1
  print(c(b, h.ASH$density[b]))
  
  #When x2 = 6.1
  x2 <- 6.1
  b <- which.min(h.ASH$breaks <= x2) - 1
  print(c(b, h.ASH$density[b]))
  
  #When x3 = 4.9
  x3 <- 4.9
  b <- which.min(h.ASH$breaks <= x3) - 1
  print(c(b, h.ASH$density[b]))
  
  #When x4 = 5.7
  x4 <- 5.7
  b <- which.min(h.ASH$breaks <= x4) - 1
  print(c(b, h.ASH$density[b]))


