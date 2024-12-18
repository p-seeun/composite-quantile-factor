library(readr) 
DATA <- read_csv("./Real data analysis/Fred-MD/2020-03.csv")
tcode = as.numeric(DATA[1,-1])
DATA = DATA[-1,] # remove row of transforming codes
rawdata = as.matrix(DATA[,-1]) #remove column of dates

# Number of series
N <- ncol(rawdata)

# Transformation function
transx = function(x, tcode){
  n = length(x)
  small = 1e-6
  y = rep(NA,n)
  
  switch(tcode,
         '1' = { 
           y <- x
         },
         '2' = { 
           y[2:n] <- x[2:n] - x[1:(n-1)]
         },
         '3' = { 
           y[3:n] <- x[3:n] - 2 * x[2:(n-1)] + x[1:(n-2)]
         },
         '4' = { 
           if (min(x, na.rm = TRUE) < small) {
             y <- NA
           } else {
             y <- log(x)
           }
         },
         '5' = { 
           if (min(x, na.rm = TRUE) > small) {
             x <- log(x)
             y[2:n] <- x[2:n] - x[1:(n-1)]
           }
         },
         '6' = { 
           if (min(x, na.rm = TRUE) > small) {
             x <- log(x)
             y[3:n] <- x[3:n] - 2 * x[2:(n-1)] + x[1:(n-2)]
           }
         },
         '7' = { 
           y1 <- rep(NA, n)
           y1[2:n] <- (x[2:n] - x[1:(n-1)]) / x[1:(n-1)]
           y[3:n] <- y1[3:n] - y1[2:(n-1)]
         }
  )
  return(y)
}

# Apply transformation to each series
data = matrix(NA, nrow = nrow(rawdata), ncol = ncol(rawdata))
for (i in 1:N) {
  data[, i] <- transx(rawdata[, i], tcode[i])
}

DATA$sasdate[c(493,732)]
data = data[493:732,]
X = data[, colSums(is.na(data)) == 0]  
X = scale(X)
dim(X) # 240 x 127

#saveRDS("./Real data analysis/Fred-MD/FRED-MD.rds")