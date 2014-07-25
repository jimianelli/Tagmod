if NOT EXIST arc mkdir arc 
echo ..\data\%1.dat > tm3.dat 
tm3 -mcmc 1000000 -mcsave 200
tm3 -mceval
copy tm3.std arc\%1.std
copy tm3.cor arc\%1.cor
copy tm3.std arc\%1.std
copy tm3.rep arc\%1.rep
copy mcout.rep arc\%1_mcout.rep
