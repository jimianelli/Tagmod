call run%1 pet_all    :: Base case
call run%1 pet_all_s1 :: Only 1 commercial and 1 charter reporting rate
call run%1 pet_all_s2 :: Free natural mortality
call run%1 pet_all_s3 :: seafisher excluded in event 8
call run%1 pet_all_s4 :: Base case with constrained M
call run%1 pet_011    :: Base case
:: call run%1 pet_011_s1 :: Only 1 commercial and 1 charter reporting rate
:: call run%1 pet_011_s2 :: Free natural mortality
:: call run%1 pet_011_s3 :: seafisher excluded in event 8
:: call run%1 pet_011_s4 :: Base case with constrained M
call run%1 seg_all
call run%1 seg_all_s1 :: constrained M
call run%1 seg_all_s2 :: M free
call run%1 seg_all_s3 :: seafisher excluded in event 8 for reporting rate
call run%1 seg_all_s4 :: Only 1 commercial and 1 charter reporting rate
call run%1 seg_011   

call run%1 tan_all
call run%1 tan_011
goto end

call run%1 tan_all_s1
call run%1 tan_all_s2
call run%1 tan_all_s3
call run%1 pet_011
call run%1 pet_011_s1
call run%1 pet_011_s2
call run%1 pet_011_s3
call run%1 seg_011
call run%1 seg_011_s1
call run%1 seg_011_s2
call run%1 seg_011_s3
call run%1 tan_011
call run%1 tan_011_s1
call run%1 tan_011_s2
call run%1 tan_011_s3
:end
