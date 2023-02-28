install.packages(dplyr)
library(dplyr)
library(hflights)
head(hflights)
flightdf <- tbl_df(hflights) #local data frame
flightdf
print(flightdf, n = 20) #this shows first 20 rows of data frame 
#this shows all
flightdf[flightdf$Month == 1 & flightdf$DayofMonth == 1] #logical expression 
#for smaller data frame, finds all flights during January
filter(flightdf, Month == 1, DayofMonth == 1) #dyplyr specific filter method
select(flightdf, Year:DayofMonth, contains("Taxi"))#dyplyr specific select method
flightdf %>% select(Year,DayofMonth) %>% filter(DayofMonth == 9) #%>% is used for 
#organzing operations better
my_db <- src_sqlite("my_db.sqlite") #connecting to SQLite database
flights_tb <- tbl(my_db, "hflights") #cnonects to hflights dataset through db
