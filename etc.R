




(pseudo-code)
row_limit = count(array);
if(row_limit > 0){
        column_limit = count(array[0]);
        for(x = max(0, i-1); x <= min(i+1, row_limit); x++){
                for(y = max(0, j-1); y <= min(j+1, column_limit); y++){
                        if(x != i || y != j){
                                print array[x][y];
                        }
                }
        }
}



x = c("d","b", "a", "c","d")
       
names(x)=  c( 1,2,3,4,5)
rank(x, ties.method="max")

x.1 <- order(unique(x))
(1,1,1,4,5)--> (1,1,1,2,3)

1.
내가 몇번째냐는 값을 얻을려면...


unique하게 해서....
rank하고.. index 값을 구한다음에...
좌표값을 가지고 비교하면서, index 값을 매핑한다...
효율적인 방법이 있을 것 같기는 하나...

unique 함수적용후 sorting한 벡터를 만들고,
해당좌표 컬럼에 match를 적용한다.(벡터연산이므로 속도가 빠름.)

1-1

rank 값을 순차적으로 올라가게 할 수 없나...
1 1 1 1 5 6 
--> 1 1 1 1 2 3 

rm(list=ls(pattern="stat1H.3"))
system.time( replicate(100,stat1H.2 <- transform(stat1H.1,id=as.numeric(factor(lon)))))
system.time( replicate(100, stat1H.3 <- transform(stat1H.1,id= match(lon, sort(unique(lon))))))


stat1H.2 <- transform(stat1H.1,xid= match(lon, sort(unique(lon))),
                               yid= match(lat, sort(unique(lat))))

xmax <- max(stat1H.2$xid)
ymax <- max(stat1H.2$yid)

row_limit = xmax;

nghborlist <-matrix(list(), nrow=xmax, ncol=ymax)
i=100; j=100;
column_limit = ymax;
for(x in seq(max(1, i-20), min(i+20, row_limit))){
        for(y in seq(max(1, j-20), min(j+20, column_limit))){
                if(x != i || y != j){
                        nghborlist[x][y] <- c(nghborlist[x][y],
                            stat1H.2[stat1H.2$xid == x & stat1H.2$yid ==y, 'cellid']);
                }
        }
}




위의 방법은 틀렸음.
x,y좌표를 인덱스화 했을뿐...x,y인덱스로 셀을 추출할때 Null 값이 많이 나옴...
처음부터 잘못 생각한 것임...

셀들의 좌표를 이용해서, 셀간 거리와 각도를 구하고, 
각도를 8개로 나눠서 그룹화하고 그중에서 거리가 최소인 셀을 인접셀로 한다. 
인접셀은 8개로 한정한다.
해당셀을 중심으로 해서 좌표의 차가 일정한 범위를 정해야 한다. 그렇지 않으면 너무 많은 셀들에 대해 계산...
- 해당 셀 좌우 5Km 정도로 먼저 해보자.. 지도상에 나타내 보면...


circleFun <- function(center = c(0,0),radius = 1.5, npoints = 100){
        # 좌표는 GPS 기준, radius는 km단위....
        
        tt <- seq(0,2*pi,length.out = npoints)
        xx <- center[1] + radius * cos(tt) / 88.9 # 1km/88.9 :GPS x (longitude) 좌표값...
        yy <- center[2] + radius * sin(tt) / 114  # 1km/114 : GPS y (latitude ) 좌표값...
        return(data.frame(x = xx, y = yy))
}

# 1.5/88.9 = 0.01687 반경 1.5 km 지역을 커버하는 원
# longitute(x좌표)의 차이는 1이 88.9Km에 해당함.
# latitute(y좌표)의 차이는 1이 114Km에 해당함.

dat <- circleFun(c(mean(cellLoc$lon), mean(cellLoc$lat)), 0.01687, npoints = 100)
#geom_path will do open circles, geom_polygon will do filled circles

p.1 + geom_path(data=dat, aes(x,y))

kor.1 <- get_map(location = c(lon = mean(cellLoc$lon), lat = mean(cellLoc$lat)),
                 maptype = "roadmap", source = "google" ,  zoom = 15)
p.1 <- ggmap(kor.1 ) 




gpsToPolarCoord <- function(x,y) {
        d <- complex(real = x*88.9, imaginary = y*114)
        r = Mod(d) 
        theta = Arg(d) %% (2*pi) * 360/(2*pi)
        list(r, theta)    # c(r, theta) 로 하면 df의 column을 입력으로 했을때 r값이 계속 나오고, 이어서 theta값이 계속 나옴...
}


findNeighbor <- function(df=stat1H.1, cellidx= "1-17124-12", radiusLimit=1.5){
        x <- df[df$cellid == cellidx,"lon"]
        y <- df[df$cellid == cellidx,"lat"]

        polarList <-gpsToPolarCoord(df$lon-x,df$lat-y)
        df$r <- polarList[[1]]
        df$theta <- polarList[[2]]
        
        df1 <- df %>% filter(r < radiusLimit, r >0) %>% 
                group_by(theta %/% 45) %>% summarize(dist.min = min(r), cellid[which.min(r)])
        c(cellid=cellidx, neighborList=df1$cellid)
}
findNeighborLists <-Vectorize(findNeighbor, vectorize.args = "cellidx")

stat1H.1.800 <- stat1H.1 %>% filter(str_detect(FREQUENCY, "800"))

neighborList.800 <-findNeighborLists(
                        stat1H.1.800, stat1H.1.800$cellid, 1.5)

length(neighborList.800)





z <- strptime("20/2/06 11:16:16.683", "%d/%m/%y %H:%M:%OS")
z # prints without fractional seconds
op <- options(digits.secs=3)
z
options(op)
z

