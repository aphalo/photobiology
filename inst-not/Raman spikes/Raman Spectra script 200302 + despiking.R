####Load more raman spectra and compare them with lierature silica phases###

library(ggplot2)
library(dplyr)
library(scales)
library(photobiology)
library(photobiologyWavebands)
library(ggspectra)
library(ggrepel)
library(tools)
library(cowplot)
library(rlang)

setwd("C:\\Users\\Agnese\\Desktop\\Reasearch\\DFG Project\\Quartz crystals\\test convert txt in csv")

### Load files and update names and header 

header <- c("Raman_shift", "Intensity")

L = list.files(".", ".txt")

NL <- file_path_sans_ext(L)

O = lapply(L, function(x) {
  spectra <- read.table(x, header = FALSE, sep = ",", dec = ".", stringsAsFactors = FALSE)
  
})

O <- lapply(O, function(x) {names(x) <- c(header) ; return(x)})

names(O) <- c(NL)

ID <- c(NL)

O <- mapply(cbind, O, "ID"=ID, SIMPLIFY = F)

### Adjust shift

O <- lapply(O, transform, Raman_shift = Raman_shift - 3.26)

### DESPIKING Whitaker and Hayes' Algorithm

  #' ###  Function to calculate modified Z Scores

  ModifiedZscore = function(x) {
  m = median(x, na.rm=TRUE)
  M = mad(x, na.rm=TRUE)
  z = (x - m) / M 
  z
  }


  #' ### Function to annihalate spikes at locations z  

  fixer = function(y, z, ma=5) {
    n = length(y) 
    yout = y
    z[1] = z[n] = 1
    spikes = which(z==1)
    for(i in spikes) {
      w = seq(max(1,i-ma),min(n,i+ma))
      w = w[z[w] == 0] 
      yout[i] = mean(y[w])
    }	
    yout 
  }
  
  #select your intensity values (y) as matrix + your raman_shift values (w) and your ID as vectors

  y <- sapply(O, function(x) x[[2]])
  w_matrix <- sapply(O, function(x) x[[1]])
  w <- w_matrix[ ,1]

  matplot(w,y, type = "l")

  #modified z-score

  z = matrix(0, nrow(y)-1, ncol(y))

  for(i in 1:ncol(y)) z[,i] = ModifiedZscore( diff(y[,i]) ) 

  z = rbind(rep(0,ncol(y)),z)

  matplot(w,z, type="l")


  #create a matrix z identifying spikes
  #' Here the user must choose a threshold, we recommend starting at a high value and decrease
  #' until optimal

  threshold = 13.6

  z = (abs(z) > threshold) * 1

  matplot(w,z, type="l")


  #' ##  Block 5: despike the spiky spectra

  spiky = seq(ncol(z)) [colSums(z) > 0] 

  for(i in spiky) y[,i] = fixer(y[,i], z[,i]) 

  matplot(w,y, type="l")

  ### add y as new column Intensity_despiked in O

  y_list <- split(y, rep(1:ncol(y), each = nrow(y)))

  O <- mapply(cbind, O, "Intensity_despiked"= y_list, SIMPLIFY = F)

###create a data frame from the list

DF_RamanSpectra <- do.call(rbind, O)

DF_RamanSpectra <- subset(DF_RamanSpectra, select=c(ID, Raman_shift, Intensity, Intensity_despiked))

###Plot data all together RAW

plot_DF_RamanSpectra <- ggplot(DF_RamanSpectra, aes(x = Raman_shift, y = Intensity, group = ID)) +
  geom_line(aes(color = ID),
            size = 1) +
  scale_color_grey() +
  scale_x_continuous(breaks = seq(0, 1500, by = 100)) +
  theme_light() +
  labs (
    x = "Raman shift (cm-1)",
    y = "Intensity")
plot_DF_RamanSpectra

###Plot data all together DESPIKING

plot_DF_RamanSpectra_despiked <- ggplot(DF_RamanSpectra, aes(x = Raman_shift, y = Intensity_despiked, group = ID)) +
  geom_line(aes(color = ID),
            size = 1) +
  scale_color_grey() +
  scale_x_continuous(breaks = seq(0, 1500, by = 100)) +
  theme_light() +
  labs (
    x = "Raman shift (cm-1)",
    y = "Intensity")
plot_DF_RamanSpectra_despiked
   
###Add info to spectra

Raman_plot_peak <- plot_DF_RamanSpectra_despiked +
                    #stat_peaks(color = "red") +
                    stat_peaks(span = 51, geom = "text", color = "red", hjust = -0.2, angle = 66, ignore_threshold = 0.1)
Raman_plot_peak

Raman_plot_peak_title <- Raman_plot_peak_despiked +
                            labs (title = "Quartz") +
                            theme(plot.title = element_text(hjust = 0.5))
Raman_plot_peak_title


#GET_PEAKS [check if you want the Intensity [2] or Intensity_despiked [4]]

O_peaks = lapply(O, function(x){
  get_peaks(x[[1]], x[[4]], 
            ignore_threshold = 0, span = 71, strict = TRUE, x_unit = "", x_digits = 3)})

header_peaks <- c("Raman_shift", "Intensity", "Label")

O_peaks <- lapply(O_peaks, function(x) {names(x) <- c(header_peaks) ; return(x)})

names(O_peaks) <- c(NL)

O_peaks <- mapply(cbind, O_peaks, "ID"=ID, SIMPLIFY = F)


