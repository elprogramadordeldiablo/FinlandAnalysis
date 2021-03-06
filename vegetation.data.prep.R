library(here)

data <- read.csv2(here("data","veg2.csv"), header=TRUE,  stringsAsFactors = T, dec = ".")

# 
 sample1 = c(23, 45, 2, 0, NA, 4)
 sample2 = c(45, 12, 3, 16, 45, 34)
# 
 newData <- list() 
 dataa = rbind(sample1, sample2)
 
 dataa
 
colnames(dataa) = c("distance", "plot", "l4m0", "l5m4", "l3m5", "l2m4")

for(i in 1:nrow(dataa)){            #
  basicInfo = dataa[i, 1:2]
  for(j in seq(3, 6, 2)){
    newSp = dataa[i, c(j, j+1)]  # newsp= par de c�lulas j (morfoesp�cie) e j+1 (altura em cm)
    if(!any(is.na(newSp)))      # se n�o for NA. 
      newData = rbind(newData,c(basicInfo, colnames(dataa)[j], newSp))
  }
}
newData

 
 
 #### test with simpler database: test_veg2.csv
 
 data <- read.csv2(here("data","test_veg2.csv"), header=TRUE,  stringsAsFactors = T, dec = ".")
 veg.table <- data.frame() 
 str(data)
 for(i in 1:nrow(data)){            #
   basicInfo = data[i, 1:70]
   for(j in seq(71, 74, 2)){
     newSp = data[i, c(j, j+1)]
     if(!any(is.na(newSp))){
       abc <- data.frame(basicInfo, colnames(data)[j], newSp)
       colnames(abc)[71:73] = " lala"
       veg.table = rbind(veg.table,abc)
     }
   }
 }
 veg.table
 
class( as.vector(c(basicInfo, colnames(data)[j], newSp)))
write.csv(veg.table, file = here("results","veg.table.csv"), row.names = TRUE)
 
class(basicInfo)
      class(colnames(data) )
            class(newSp)
 
            #### test with the real database
            
            real.data <- read.csv2(here("data","veg2.csv"), header=TRUE,  stringsAsFactors = T, dec = ".")
            str(real.data)
            veg.table.real <- data.frame() 
          
            for(i in 1:nrow(real.data)){            #
              basicInfo2 = real.data[i, 1:71]
              for(j in seq(72, 870, 2)){
                newSp.real = real.data[i, c(j, j+1)]
                if(!any(is.na(newSp.real))){
                  abcd <- data.frame(basicInfo2, colnames(real.data)[j], newSp.real)
                  colnames(abcd)[72:74] = " lala"
                  veg.table.real = rbind(veg.table.real,abcd)
                }
              }
            }
            veg.table.real
            
            write.csv(veg.table.real, file = here("results","veg.table.real.csv"), row.names = TRUE)
            
        
 
 
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






