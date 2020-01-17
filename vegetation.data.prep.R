library(here)

data <- read.csv2(here("data","veg2.csv"), header=TRUE,  stringsAsFactors = T, dec = ".")

# 
 sample1 = c(23, 45, 2, 0, NA, NA)
 sample2 = c(45, 12, 3, 16, 45, 34)
# 
 newData <- list() 
 data = rbind(sample1, sample2)
# colnames(data) = c("distance", "plot", "l4m0", "l5m4", "l3m5", "l2m4")
# 
# for(i in 1:nrow(data)){            # 
#   basicInfo = data[i, 1:2]
#   for(j in seq(3, 6, 2)){
#     newSp = data[i, c(j, j+1)]
#     if(!any(is.na(newSp)))
#       newData = rbind(newData,c(basicInfo, colnames(data)[j], newSp))
#   }
# }
# newData

 
 
 #### test with simpler database: test_veg2.csv
 
 data <- read.csv2(here("data","test_veg2.csv"), header=TRUE,  stringsAsFactors = T, dec = ".")
 veg.table <- data.frame() 
 
 for(i in 1:nrow(data)){            #
   basicInfo = data[i, 1:70]
   for(j in seq(71, 870, 2)){
     newSp = data[i, c(j, j+1)]
     if(!any(is.na(newSp)))
       veg.table = rbind(veg.table,c(basicInfo, colnames(data)[j], newSp))
   }
 }
 veg.table
 
 
 
 
 
 
 
 
 
 
 
#### original script backup
data <- read.csv2(here("data","veg2.csv"), header=TRUE,  stringsAsFactors = T, dec = ".")
veg.table <- data.frame() 

for(i in 1:nrow(data)){            #
  basicInfo = data[i, 1:70]
  for(j in seq(71, 870, 2)){
    newSp = data[i, c(j, j+1)]
    if(!any(is.na(newSp)))
      veg.table = rbind(veg.table,c(basicInfo, colnames(data)[j], newSp))
  }
}
veg.table






