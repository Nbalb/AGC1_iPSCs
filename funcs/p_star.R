# Function that takes in a vector of p-values and outputs a corresponding number of stars 
p_star <- function(x){
  
  if(sum(is.na(x)) != 0){
    
    warning("Vector contains NA values")
  }
  
  star <- factor()
  star <- sapply(x, function(i) {
                   if(i <= 10^(-4)){"****"
                   }else if(i <= 10^(-3)){"***"
                   }else if(i <= 10^(-2)){"**"
                   }else if(i <= 5*10^(-2)){"*"
                   }else{""}
                 })

  return(star)
  
}
