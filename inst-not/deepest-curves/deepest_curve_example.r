# Functional data analysis of multiple spectra, i.e. "finding the deepest curve"
library(fda)
library(fda.usc)

Beo1leaf00000_irrad <- read.table("inst-not/deepest-curves/Beo1leaf00000_irrad.txt", header = TRUE)
# Beo1leaf00000_irrad is a Maya scan from Lammi, 22.5.2015, from a Betula forest
  # There are 100 scans under leaf shade, i.e. under the Betula canopy
# in a spot with shade from the leaves above
# First the data needs to be transposed
Beo1leaf00000_irrad_t <-t(Beo1leaf00000_irrad[,-1])

fdataobj <- fdata(Beo1leaf00000_irrad_t)

# Returns the deepest curve following FM criteria
func_med_FM <- func.med.FM(fdataobj)
# Returns the deepest curve following mode criteria
func_med_mode <- func.med.mode(fdataobj)
# Returns the deepest curve following RP criteria
func_med_RP <- func.med.RP(fdataobj)

# There are further options, please see the fda.usc package and
# the article about the package
# https://www.jstatsoft.org/article/view/v051i04/v51i04.pdf

# Here is a plot to see where the deepest curve really is
# Comparing both methods above
graphics.off()
# Guessing the ymax limit here
plot(func_med_FM, ylim=c(0,0.3), main=NULL)
# The loop plots the original data
for (i in 1:nrow(Beo1leaf00000_irrad_t))
{
  lines(Beo1leaf00000_irrad_t[i,])
}
# Then the deep lines are plotted in different colours to see "which depth to use"
# med_FM and med_RP give so similar results here that the lines overlap
lines(func_med_FM, col=3) # median FM = green
lines(func_med_mode, col=4) # median mode = blue
lines(func_med_RP, col=2) # median RP = red

# Saving in .txt format only the deepest curve
irrad.out <- Beo1leaf00000_irrad
deepest.irrad.out.file.name <- paste(irrad.out.file.name,"_deepestirrad.txt",sep="")
write.table(irrad.out[,c(1,fda.out$medcurve+1)], deepest.irrad.out.file.name, row.names=FALSE)


# The following can be tried if e.g. for part of the scans the sensor
# has been saturated, as it has here, in a measurement done in a sunfleck
# It applies the outliers detection method based on trimming with mode depth
# Please see ?depth.FM and ?depth.mode and ?depth.RP for explanation of the difference between them
# For the example, suggested to use depth.FM, just based on looking the graph "which depth to use"

Beo1sun00000_irrad_t <-t(Beo1sun00000_irrad[,-1])

fdataobj <- fdata(Beo1sun00000_irrad_t)

out.trim.FM <-  outliers.depth.trim(fdataobj, dfunc = depth.FM, nb = 20, smo = 0.1, trim = 0.06)
func_med_FM <- func.med.FM(fdataobj[out.trim.FM[[1]]])

out.trim.mode <-  outliers.depth.trim(fdataobj, dfunc = depth.mode, nb = 20, smo = 0.1, trim = 0.06)
func_med_mode <- func.med.mode(fdataobj[out.trim.mode[[1]]])

out.trim.RP <-  outliers.depth.trim(fdataobj, dfunc = depth.RP, nb = 20, smo = 0.1, trim = 0.06)
func_med_RP <- func.med.RP(fdataobj[out.trim.RP[[1]]])

# Graph "which depth to use"
plot(fdataobj, col=1) # black
lines(func_med_FM, col=3) # green
lines(func_med_mode, col=4) # blue
lines(func_med_RP, col=2) # red
