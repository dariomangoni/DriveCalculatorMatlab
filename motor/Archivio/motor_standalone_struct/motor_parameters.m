function [deviazione, R, L] = motor_parameters(U, omega, I, Ke, n0)
    Z = (U - Ke*omega)./I;
    
    bench_len = size(Z,1);
    
    R_ind = 0;
    for m = 1:bench_len-1
        for n=m+1:bench_len
            if (~n0)&&((m==1)||(n==1))
                continue
            end
            R_ind = R_ind+1;
            R_temp(R_ind) = sqrt(abs((Z(n)^2*omega(m)^2-Z(m)^2*omega(n)^2 )/(omega(m)^2-omega(n)^2)));
            L_temp(R_ind) = sqrt(abs((Z(n)^2-Z(m)^2)/(omega(m)^2-omega(n)^2)));
        end
    end
    
    R = mean(R_temp);
    L = mean(L_temp);
    I_calc = (U - Ke*omega) ./ sqrt(R^2 + (L*omega).^2);
    
    deviazione = sum(sqrt((I_calc-I).^2./(R_ind*I.^2)),1);


end