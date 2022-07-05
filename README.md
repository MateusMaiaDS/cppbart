# cppbart
Implementation in C++ in R


``` r
# install.packages("devtools")
devtools::install_github("MateusMaiaDS/cppbart")
```


### Simple example

```r 
ibrary(gprcpp)
n <- 500
n_new <- 200
# Creating a simple simulation scenario
x <- matrix(seq(-pi,pi,length.out = n))
x_new <- matrix(seq(-pi,pi,length.out = n_new))
y <- sin(x) + rnorm(n = n,sd = 0.5)

gp_mod <- gp_main(x = x,x_star = x_new,y = y,
                  gp_param = list(phi=0.1,nugget = 0.25, nu = 1),phi_sample = T)


# Plotting the GP
plot(x,y,pch=20)
lines(x_new, gp_mod$y_hat, col = "orange")
lines(x_new, gp_mod$y_hat + 1.96*gp_mod$y_sd, lty = "dashed" , col  = "orange")
lines(x_new, gp_mod$y_hat - 1.96*gp_mod$y_sd, lty = "dashed" , col = "orange")
lines(x_new, sin(x_new), col = "blue")

```
