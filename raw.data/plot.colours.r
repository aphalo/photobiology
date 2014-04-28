colCount <- 25 # number per row
rowCount <- trunc(length(my.colors) / colCount) + 1

plot( c(1,colCount), c(0,rowCount), type="n", ylab="", xlab="",
      axes=FALSE, ylim=c(rowCount,0))
title("400 to 700 nm colors")

for (j in 0:(rowCount-1))
{
  base <- j*colCount
  remaining <- length(my.colors) - base
  RowSize <- ifelse(remaining < colCount, remaining, colCount)
  rect((1:RowSize)-0.5,j-0.5, (1:RowSize)+0.5,j+0.5,
       border="black",
       col=my.colors[base + (1:RowSize)])
#  text((1:RowSize), j, paste(base + (1:RowSize)), cex=0.7,
#       col=TextContrastColor[base + (1:RowSize)])
}
