set path=c:\Rtools\bin;c:\Rtools\perl\bin;c:\Rtools\gcc-4.6.3\bin;c:\R\bin;c:\Program Files\R\R-3.0.1\bin\i386
cd C:/Users/%USERNAME%/Documents/Rpackages
Rcmd.exe INSTALL --build --library=C:/Users/%USERNAME%/Documents/R/win-library/3.0 photobiology_0.2.8.tar.gz
cd photobiology