## COVERAGE PLOTS USING PLOTLY AND BASE R
## Author :  Mihai Croicu, Uppsala University, 2018.
## License : https://opensource.org/licenses/MIT
## For AQM 2019


##########################
#PREPPING
##########################


#Clean up workspace of old objects
rm(list=ls())

#Install plotly and processx if you don't have it.
if (!require("processx")) install.packages("processx")
if (!require("plotly")) install.packages("plotly")
if (!require("webshot")) {
  install.packages('webshot')
  webshot::install_phantomjs()
  }

#Set seed to 0 for replication
set.seed(0)

#We're going to use plotly for countour plots today.
library(plotly)


##########################
#SETTING UP THE DGP AND SAMPLING Xs
##########################

#We are going to plot the real DGP, not a simulated model.
#The formula below is the one we use:
#If you want to plot simulation results you are going to need to extract the coefficients
#The formula we use is
#   y_hat = 0.2 + 0.5*x1 + 0.25*x2 - 0.1*x1*x2
b0 <- .2 
b1 <- .5 
b2 <- .25
b12 <- -.1

# We are going to draw n times X1s and X2s.
# X1 and X2 are random variables drawn from the normal distribution.
# With mean=3 and sd=1
# For the contour plot, the algorithm will need to calculate n^2 values!
# As it needs to get the y_hat at each combination of X1 and X2 sampled!

n <- 500
X1 <- rnorm(n, 3, 1)
X2 <- rnorm(n, 3, 1)


##########################
#PREPARING DATA FOR THE COVERAGE/3D PLOT
##########################


# We need X1 and X2 sorted for plotting purposes.
# It doesn't matter that we break the connection between X1 and X2 in the data, 
# as all combinations of X1 and X2 are tested (all n * n combinations),
# including those in the paired data.

X1_sort <- sort(X1)
X2_sort <- sort(X2)

# We have to write a function that takes ONE x1 value (x1_i), ONE x2 value (x2_i) and produces ONE y_hat (y_hat_i) value
# We will then run this function for all possible combinations of X1_sort and X2_sort to produce all possible interaction values.
# That's wy we use x1 and x2 with small x, so that there's no confusion to the X1 vectors!

z.func <- function(x1,x2) {b0 + b1*x1 + b2*x2 + b12*x1*x2}

# Try it! 
# z.func(2,3) yields 1.35
# z.func(5,1.5) yields 2.325

# We now create, for each X1 and X2 a matrix of Y_hats.
# We apply the formula from X1 with each element of X2, storing each multiplication of X1_i as its own rows
# Like this:
# X1 = [X1_1, X1_2, X1_3,...., X1_N]
# X2 = [X2_1, X2_2, X2_3,...., X2_N]
#
#          COL_1         COL_2        COL_3           COL_4        .. .. ..   COL_N
# ROW_1  f(X1_1,X2_1) f(X1_1,X2_2)  f(X1_1,X2_3)   f(X1_1,X2_4)              f(X1_1,X2_N)
# ROW_1  f(X1_2,X2_1) f(X2_1,X2_2)  f(X1_2,X2_3)   f(X1_2,X2_4)              f(X1_2,X2_N)
# ROW_1  f(X1_3,X2_1) f(X3_1,X2_2)  f(X1_3,X2_3)   f(X1_3,X2_4)              f(X1_3,X2_N)
# ...        :
# ...        :              
# ROW_N  f(X1_N,X2_1) f(X1_N,X2_2) f(X1_N,X2_3)    f(X1_N,X2_4)              f(X1_N,X2_N)
#
# This is called the outer product of two vectors
# And also known as the matrix multiplication of a 1N to a N1 matrix 
# And also known as the tensor multiplication of two 1D tensors.
# And R knows how to do this by default, because matrix multiplication is essential for many things.
# This is done using the outer() function.
# This takes the two vectors as the first two inputs, and the function working on individual elements X1_i and X2_i
# That's why we needed to define the functions above in terms of elements, not vectors!

z.matrix <- outer(X1_sort,X2_sort, z.func)

### EXTRABITS!!!
### We can also do this by hand!
### Here's how:
### 1. Preallocate memory for a n*n matrix. We fill it with the index, just because it's fastest
#z.matrix_by_hand = matrix(1:(n*n),nrow=n)
### i indexes rows and j columns, thus z.matrix_by_hand[i,j] is row i, column j of the finished matrix.
#for (i in 1:n) {
#  for (j in 1:n) {
#    z.matrix_by_hand[i,j] = z.func(X1_sort[i],X2_sort[j])
#  }
#}
### We can check whether the two matrixes (the one made by hand and by outer are the same)
# all(z.matrix_by_hand==z.matrix) # TRUE
### END OF EXTRABITS!!!

##########################
#PLOTTING WITH PLOTLY
##########################

# Plotly wants to put X1 and X2 "flipped" from what we expect (X1 on the vertical and X2 on the horizontal)
# We instead transpose the z.matrix to get the coordinates right. 
# This is only a plotly thing, If you use other plotters, you won't need the transposition.
transposed.z.matrix <- t(z.matrix)

## And here we plot.
contour = plot_ly(x=X1_sort, y=X2_sort, z=transposed.z.matrix, type='contour') %>% layout(xaxis = list(title='X1'), yaxis = list(title='X2'))
surface = plot_ly(x=X1_sort, y=X2_sort, z=transposed.z.matrix) %>% add_surface()

## THESE PLOTS ARE INTERACTIVE! TRY TO SPIN THE SURFACE PLOT!
## Cool, eh?
contour
surface

## And we can also save these as pdfs. You will need to install the plotly orca server or go through their server.
## Don't worry, if you choose to do these with Plotly, just leave these commented, we can check on our own machine.
## If you need help exporting a plotly plot to your machine, for example, for your own work, contact us!
export(contour,'contour.pdf')
export(surface,'surface.pdf')


##########################
#PLOTTING WITH BASE R
##########################

## Uncomment below code to work with base-R. 
## These plots are quite ugly.
## Observe that you don't need the transposition for the Z matrix when using base plots

##pdf("contour_plot_base.pdf", width = 64, height = 64) 
##cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(24))
##filled.contour(x=X1_sort, y=X2_sort, z=z.matrix, xlab="X1", ylab="X2", 
##               col = cols)
##dev.off() 

##pdf("persp_base.pdf", width = 64, height = 64) 
##persp(x=X1_sort, y=X2_sort, z=z.matrix)
##dev.off() 

