if NOT EXIST arc mkdir arc 
echo ..\data\%1.dat > tm3.dat 
tm3 -nox
copy tm3.std arc\%1.std
copy tm3.cor arc\%1.cor
copy tm3.std arc\%1.std
copy tm3.rep arc\%1.rep
:: tm3 -ind 2011only_od.dat -mceval
:: tm3 -ind alldat_od.dat -mcmc 1000000 -mcsave 200
:: tm3 -ind alldat_od.dat -mceval
:: copy tm3.psv alldat.psv
cleanad
