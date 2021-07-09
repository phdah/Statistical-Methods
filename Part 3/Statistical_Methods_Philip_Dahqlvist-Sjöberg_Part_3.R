# Title     : Statistical Method Part 3
# Created by: Philip Dahqlvist-SJöberg
# Created on: 20/06/02
table_value <- matrix(c(0, 	0	,92000 ,	20000 ,	0
,22000 ,	1000 ,	0 ,	0 ,	0,
1238000 ,	58000 ,	97000 ,	220000 ,	0,
63000 ,	146000 ,	112000 ,	495000 ,	0
), ncol = 5, byrow = T)

table_largest <- matrix(c(0, 	0	,53000, 	20000, 	0,
16000 ,	1000 ,	0, 	0, 	0,
178000 ,	23000 ,	57000 ,	35000 ,	0,
21000 ,	71000 ,	87000 ,	134000 ,	0
), ncol = 5, byrow = T)

table_sec_largest <- matrix(c(0, 	0,	34000 ,	0, 	0,
6000 ,	0,0,0,0,
74000 ,	18000 ,	39000 ,	25000 ,	0,
15000 ,	69000 	,14000 ,	132000 ,	0
), ncol = 5, byrow = T)

# 2a, dominance (1,60)
dominance1 <- function (table_largest, table_value, percent){
  output <- round(table_largest / table_value, 4)
  output[is.nan(output)] <- 0
  output[output > percent/100] <- "sdc"
  return(output)
}
output <- dominance1(table_largest, table_value, 60)
output

# 2a, dominance (2,90)
dominance2 <- function (table_largest, table_sec_largest, table_value, percent){
  output <- round((table_sec_largest + table_largest) / table_value, 4)
  output[is.nan(output)] <- 0
  output[output > percent/100] <- "sdc"
  return(output)
}
output <- dominance2(table_largest, table_sec_largest, table_value, 90)
output

# 2b, p%-rule (1,11)
p_percent <- function (table_largest, table_sec_largest, table_value, percent){
  output <- table_value - table_sec_largest
  output_tmp <- output - table_largest
  output <- round(output_tmp / output, 4)
  output[output < percent/100] <- "sdc"
  output[output == "NaN"] <- 0
  return(output)
}
output <- p_percent(table_largest, table_sec_largest, table_value, 11)
output