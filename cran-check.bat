setlocal
set PATH=C:\Program Files\R\R-devel\bin\x64;%PATH%
R CMD check --as-cran ../photobiology_0.9.4.tar.gz
pause
