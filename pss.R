library(readxl)
library(ggplot2)
library(stringr)
library(dplyr)

glrvs <- guides(fill=guide_legend(reverse = TRUE), color=guide_legend(reverse = TRUE)) 

setwd("~/Project_timeSeriesAnalysis/data")
x=list.files(pattern = "DLU.*")
x
pssData = read_excel("DLUL_PRB_TOT_PRB_1110537_20160919105501.xls",skip=2)
dim(pssData)
colnames(pssData)[1] <- "complexinfo"

# strdf <- pssData[,1]  # dataframe

pssData$eNB <- str_extract(pssData$complexinfo, "\\d-\\d{5}")
#pssData$eNB <- str_extract(pssData[,1], "\\d-\\d{5}")  #cannot manage dataframes as inputs..
selDeleteLine <- str_detect(pssData$complexinfo, "[소총] 계")
pssData <- pssData[!selDeleteLine,]
dim(pssData)



colnames(pssData) <- str_replace(colnames(pssData),"\\(.*\\)", "")
ggplot(pssData, aes(x=TotPrbDLMax, fill= FREQUENCY)) + geom_histogram(binwidth = 1) + glrvs

heavyCell  <- pssData[pssData$TotPrbDLMax >90,]
ggplot(heavyCell, aes(x=TotPrbDLMax, fill= FREQUENCY)) + geom_histogram(binwidth = 1)
dim(heavyCell)
heavyCellList <- pssData %>% select(eNB, CNUM)


#################
list.files()



pss5min <- read.csv("pss5min1.csv")
selDeleteLine = str_detect(pss5min[,2], "[소총] 계")
pss5min <- pss5min[!selDeleteLine, ]

pss5min$dist <- cut(pss5min$TotPrbDLMax,breaks = seq(0,100, by=10),labels = seq(10,100, by=10), include.lowest = TRUE)



ggplot(pss5min, aes(x=TotPrbDLMax, fill= FREQUENCY)) + geom_histogram(binwidth = 10) + glrvs
ggplot(pss5min, aes(x=dist, fill= FREQUENCY)) + geom_bar() + glrvs


summary(pss5min$dist)
plot(pss5min$dist)

table(pss5min$FREQUENCY)
pss5min  %>% filter(TotPrbDLMax < 1.) %>% select(FREQUENCY) %>% table
# 2.6GHz not activated..
unique(pss5min$FREQUENCY)
selDelLine <- str_detect(pss5min$FREQUENCY, "2.6GHz")
pss5min <- pss5min[!selDelLine, ]


# to find cell id on which max prb is above 90
heavyList5min <- pss5min %>%  filter(TotPrbDLMax > 90) %>% select(eNB, CNUM, time,TotPrbDLMax)

# time series graph...
ggplot(pss5min %>% filter(eNB %in% c('1-26507','1-26508')), aes(x=time, y=TotPrbDLMax, group=CNUM)) + geom_line()
ggplot(pss5min %>% filter(eNB =='1-26507'), aes(x=time, y=TotPrbDLMax, color=CNUM, group=CNUM)) + geom_line()
ggplot(pss5min %>% filter(eNB =='1-26507', FREQUENCY=="800MHz"), aes(x=time, y=TotPrbDLMax, color=CNUM, group=CNUM)) + geom_line() + ggtitle("800MHz")
ggplot(pss5min %>% filter(eNB =='1-26507', FREQUENCY=="1.8GHz(20M)"), aes(x=time, y=TotPrbDLMax, color=CNUM, group=CNUM)) + geom_line()+ ggtitle("1.8GHz")
ggplot(pss5min %>% filter(eNB =='1-26507', FREQUENCY=="2.1GHz(20M)"), aes(x=time, y=TotPrbDLMax, color=CNUM, group=CNUM)) + geom_line()+ ggtitle("2.1GHz")

##########################################

mungfile <- function(filename,date){ # date:string, "20160920"
        var_df <- read.csv(filename)
        # selDeleteLine = str_detect(var_df[,2], "[소총] 계")
        # var_df <- var_df[!selDeleteLine, ]
        var_df$dist <- cut(var_df$TotPrbDLMax,breaks = seq(0,100, by=10),
                           labels = seq(10,100, by=10), include.lowest = TRUE)
        var_df$datetime <- str_c(date, var_df$time, sep=" ")
        var_df$datetime <- as.POSIXct(strptime(var_df$datetime,"%Y%m%d %H:%M", tz="GMT"))
        var_df$date <- date
        cellnum <- str_extract(var_df$CNUM, "\\d\\d?")
        var_df$cellid <- str_c(var_df$eNB, cellnum, sep = "-")
        return(var_df)
}

lfiles <- list.files(pattern = "pss5M\\d{8}")
for(file in lfiles){
        fdate <- str_extract(file, "\\d{8}")
        vdate <- str_sub(fdate,5)
        assign(paste0("p5M", vdate), mungfile(file,fdate))
}
##########################

# 
# str_extract(lfiles[1], "\\d{4}(\\d{4})")
# fff <- str_extract(lfiles[1], "\\d{4}(\\d{4})")
#         mungfile(file,fdate)
# p5min0919 <- mungfile("pss5min0919.csv","2016-09-19")
# p5min0920 <- mungfile("pss5min0920.csv","2016-09-20")

p5min <-rbind(p5min0919, p5min0920)

ggplot(p5min0919, aes(x=dist, fill= FREQUENCY)) + geom_bar() + glrvs
ggplot(p5min0920, aes(x=dist, fill= FREQUENCY)) + geom_bar() + glrvs



ggplot(p5min,aes(x=dist, fill= FREQUENCY)) + geom_bar() + glrvs


ggplot(p5min %>% filter(eNB =='1-26511', FREQUENCY=="800MHz"), aes(x=time, y=TotPrbDLMax, color=CNUM, group=CNUM)) + geom_line(size=2) + facet_wrap(~date, ncol = 1)

#########
# to get some matric but not useful...

library(dtw)
ref = p5min %>% filter(eNB =='1-26511', FREQUENCY=="800MHz", CNUM=='cNum0', date=="2016-09-19") %>% select(TotPrbDLMax)
ref <- ref[ , 1]
query = p5min %>% filter(eNB =='1-26511', FREQUENCY=="800MHz", CNUM=='cNum0', date=="2016-09-20") %>% select(TotPrbDLMax)
query <- query[ , 1]
alignment <- dtw(query, ref, step = symmetric2, keep = TRUE, open.end = TRUE)
alignment$normalizedDistance
alignment$costMatrix
alignment$directionMatrix
alignment$stepPattern
alignment$index1
alignment$index2
alignment$localCostMatrix
alignment$distance
#######################################
alignment <- dtw(query, ref, step = symmetric1)
alignment$distance
alignment$index1
alignment$index2
alignment$stepsTaken
alignment$stepPattern
#######################################

mean(ref)
var(ref)
mean(query)
var(query)

#######################################

# get  time series of each cell on a matrix..

mkMTX <- function(df, Freq){
        cellList <- df %>% filter(FREQUENCY==Freq) %>% select(cellid)
        cellList <- cellList[ ,1] # becomes a vector.
      
        mtx <- matrix(data=NA,nrow=length(cellList),ncol=12)
        for(i in 1:length(cellList)){
                selCellPrbMax <- 
                        df %>% filter(cellid==cellList[i], FREQUENCY==Freq) %>% 
                        select(date, time, TotPrbDLMax)
                mtx[i,] <- selCellPrbMax$TotPrbDLMax
                
        }
        return(mtx)
}

mtx0919 <- mkMTX(p5M0919, "800MHz")
mtx0920 <- mkMTX(p5M0920, "800MHz")
mtx0921 <- mkMTX(p5M0921, "800MHz")

diff_mtx <- mtx0919 - mtx0920

m19mean <- apply(mtx0919,1,mean)
m19max <- apply(mtx0919,1,max)

hist(m19max)
hist(m19mean)

m20mean <- apply(mtx0920,1,mean)
m20max <- apply(mtx0920,1,max)

hist(m20max)
hist(m20mean)

diff_mean <- apply(abs(diff_mtx),1,mean)
diff_var <- apply(abs(diff_mtx),1,var)

hist(diff_mean)
hist(diff_var)

arr_prb <- array(0, c(length(cellList), 12, 5))

#
arr_prb[,,1] <- mtx0919
arr_prb[,,2] <- mtx0920
arr_prb[,,3] <- mtx0921


