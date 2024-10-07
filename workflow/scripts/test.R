# Load necessary library
library(ggplot2)

# Simulate some data
set.seed(123) # for reproducibility
data_sim <- data.frame(
  var1 = rnorm(100),  # 100 random normal values for var1
  var2 = rnorm(100)  # 100 random normal values for var2
)
#${userHome}/micromamba/envs/r/bin/R
#${userHome}/micromamba/envs/r/bin/radian
# View first few rows of the dataset
head(data_sim)

# Plot the data using ggplot2
ggplot(data_sim, aes(x = var1, y = var2)) +
  geom_point() +  # create a scatterplot
  theme_minimal() +  # use minimal theme for the plot
  ggtitle("Scatterplot of Simulated Data") +  # add a title to the plot
  xlab("Variable 1") +  # label for x-axis
  ylab("Variable 2")  # label for y-axis


# Defining the Fibonacci function
fibonacci <- function(n) {
  if(n <= 1){
    return(n)
  } else {
    a <- fibonacci(n-1) + fibonacci(n-2)
    return(a)
  }
}

# Calculating and printing fibonacci sequence and the computation time required
start_time <- Sys.time()  # Get the system time just prior to starting the calculations

# Calculating the 35th term in the Fibonacci sequence
print(fibonacci(10))

end_time <- Sys.time()  # Get the system time immediately after the calculations are done

# Printing the calculation time
print(end_time - start_time)