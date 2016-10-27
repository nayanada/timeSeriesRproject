interleave <- function(v1,v2)
{
        ord1 <- 2*(1:length(v1))-1
        ord2 <- 2*(1:length(v2))
        c(v1,v2)[order(c(ord1,ord2))]
}

interleave(4:13, 14:23)

# 1. 원하는 순서를 해당 데이터에 맵핑시키고...
# 2. 적은 순서에 해당되는 것을 order 함수에 적용하면.... 원하는 순서대로 해당 인덱스 값을 배열해준다.
# 3. 인덱스 값을 배열에 다시 적용하면 된다....


(x <- c(sort(sample(1:20, 9)), NA))
(y <- c(sort(sample(3:23, 7)), NA))
union(x, y)
intersect(x, y)
setdiff(x, y)
setdiff(y, x)
setequal(x, y)

