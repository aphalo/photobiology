library(photobiology)

wl <- c(350, 800)

# my.colors <- w_length2rgb(wl[1]:wl[2])
# my.colors <- w_length2rgb(wl[1]:wl[2], sens=ciexyzCC2.data)
my.colors <- w_length2rgb(wl[1]:wl[2], sens=ciexyzCMF2.data)

colCount <- 25 # number per row
rowCount <- trunc(length(my.colors) / colCount) + 1

plot( c(1,colCount), c(0,rowCount), type="n", ylab="", xlab="",
      axes=FALSE, ylim=c(rowCount,0))
title(paste("RGB colours for", as.character(wl[1]), "to", as.character(wl[2]), "nm"))

for (j in 0:(rowCount-1))
{
  base <- j*colCount
  remaining <- length(my.colors) - base
  RowSize <- ifelse(remaining < colCount, remaining, colCount)
  rect((1:RowSize)-0.5,j-0.5, (1:RowSize)+0.5,j+0.5,
       border="black",
       col=my.colors[base + (1:RowSize)])
}

w_length2rgb(580)
w_length2rgb(580, sens=ciexyzCC.data)
w_length2rgb(580, sens=ciexyz2006CMF.data)

w_length2rgb(395)
w_length2rgb(395, sens=ciexyzCC.data)
w_length2rgb(395, sens=ciexyz2006CMF.data)

w_length2rgb(750)
w_length2rgb(750, sens=ciexyzCC.data)
w_length2rgb(750, sens=ciexyz2006CMF.data)
