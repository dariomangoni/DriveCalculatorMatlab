function [U, I, eta] = motor_function( omega, Mout, motor_id )

% Required data - from propeller
% omega = 9500*pi/30;  % required angular velocity
% Mout = 0.11; % required output torque

% Setted by optimizator or by user
% motor_id = 27;
prefer_benchmark_over_rough_estimate = true; % TRUE if we want to extract data from database instead of using rough estimation*
extract_ALWAYS_parameters_from_benchmark = false; % TRUE if we always want to extract data from database
%* rough estimation is when we have only no load data and some parameters
%without any on-load benchmark

% Benchmark extraction from database
load('database.mat')
mymotor_row = find(motor_table(:,1)==motor_id); %#ok<*NODEF>
if (isempty(mymotor_row)) || (length(mymotor_row)>1)
    disp(['Motor ID ', num2str(motor_id),' ambiguous']);
    U = NaN;
    I = NaN;
    eta = NaN;
    return
end

bench_table = mdata_table(mdata_table(:,3)==motor_id,:); % separates only the rows of "mdata_table" under analysis
U0 = bench_table(bench_table(:,7)==1,4);
omega0 = bench_table(bench_table(:,7)==1,5)*pi/30;
I0 = bench_table(bench_table(:,7)==1,6);
Ub = bench_table(bench_table(:,7)~=1,4);
omegab = bench_table(bench_table(:,7)~=1,5)*pi/30;
Ib = bench_table(bench_table(:,7)~=1,6);

if extract_ALWAYS_parameters_from_benchmark || ((size(Ub,1)==0)&&prefer_benchmark_over_rough_estimate)% if we have only one noload data set
    kV = motor_table(mymotor_row,15); % rpm/V
    R = motor_table(mymotor_row,16)*1.66; %ohm
    Ke = 1/(kV*pi/30); %Nm/A or V/(rad/s)
    L = motor_table(mymotor_row,17); % ohm/(rad/s)
    MvK = motor_table(mymotor_row,18); % Nm
    kL = motor_table(mymotor_row,19); % Nm/(rad/s)
elseif (size(Ub,1)==0)
        R = motor_table(mymotor_row,16)*1.66; %ohm
        nsid = (U0(1) * kV)/(U0(1)-R*I0(1));
        Ke = 30/pi/nsid;
        L = 0; % ohm/(rad/s)
        MvK = 0; % Nm
        kL = I0(1)*Ke/omega0(1); % Nm/(rad/s)
else % there is almost one other data set
    
    n0=true;
    test = 0;
    Ke_mean = 0;
    
    % check if there are omegas for at least one noload data and create a
    % temporary benchmark data array
    if (omega0(1) == 0)
        if exist('omega0(2)','var')&&(omega0(2)~=0)
            U_temp = [U0(2); Ub]; omega_temp = [omega0(2); omegab]; I_temp = [I0(2); Ib];
        else
            n0 = false;
            U_temp = [U0(1); Ub]; omega_temp = [omega0(1); omegab]; I_temp = [I0(1); Ib];
        end
    else
    U_temp = [U0(1); Ub]; omega_temp = [omega0(1); omegab]; I_temp = [I0(1); Ib];
    end
    bench_len = size(U_temp,1);
       
    
    % test with all data
    if (n0==true)
        test = test+1;
        Ke_arr(test) = (sum(I_temp.*U_temp,1)*sum(omega_temp.*I_temp,1) - sum(omega_temp.*U_temp,1)*sum(I_temp.*I_temp,1))/...
                 (sum(omega_temp.*I_temp,1)*sum(omega_temp.*I_temp,1) - sum(omega_temp.*omega_temp,1)*sum(I_temp.*I_temp,1));
        R_arr(test) = (sum(I_temp.*U_temp,1)*sum(omega_temp.*omega_temp,1) - sum(omega_temp.*U_temp,1)*sum(omega_temp.*I_temp,1))/...
                 (sum(I_temp.*I_temp,1)*sum(omega_temp.*omega_temp,1) - sum(omega_temp.*I_temp,1)*sum(omega_temp.*I_temp,1));
        L_arr(test) = 0;

        I_calc = (U_temp - Ke_arr(test)*omega_temp) ./ R_arr(test);

        deviazione(test) = sum(sqrt((I_calc-I_temp).^2./(bench_len.*I_temp.^2)),1);
    end
    
    % test with two set of data
    N_test = 0;
    for m = 1:bench_len-1
        for n=m+1:bench_len
            if (~n0)&&((m==1)||(n==1))
                continue
            end
            test = test+1;
            N_test =N_test+1;
            Ke_arr(test) = (U_temp(m)*I_temp(n) - U_temp(n)*I_temp(m)) / (I_temp(n)*omega_temp(m) - I_temp(m)*omega_temp(n));
            [deviazione(test), R_arr(test), L_arr(test)] = motor_parameters(U_temp, omega_temp, I_temp, Ke_arr(test), n0);
            Ke_mean = Ke_mean + Ke_arr(test);
        end
    end
    
    % test with mean of two data set
    Ke_mean = Ke_mean/N_test;
    test = test+1;
    Ke_arr(test) = Ke_mean;
    [deviazione(test), R_arr(test), L_arr(test)] = motor_parameters(U_temp, omega_temp, I_temp, Ke_arr(test), n0);
    
    % test with three data set
    if (bench_len>=4)
        Ke_mean = 0;
        N_test = 0;
        for m = 1:2
            for n = m+1:3
                for p = n+1:4
                    if (~n0)&&((m==1)||(n==1)||(p==1))
                        continue
                    end
                    test = test+1;
                    Ke_arr(test) = - (- I_temp(m) * I_temp(p) * U_temp(n) * omega_temp(m)  +  I_temp(m) * I_temp(n) * U_temp(p) * omega_temp(m)  +  I_temp(n) * I_temp(p) * U_temp(m) * omega_temp(n)...
                                    -  I_temp(m) * I_temp(n) * U_temp(p) * omega_temp(n)  -  I_temp(n) * I_temp(p) * U_temp(m) * omega_temp(p)  +  I_temp(m) * I_temp(p) * U_temp(n) * omega_temp(p)) /...
                                    (I_temp(m) * I_temp(p) * omega_temp(m) * omega_temp(n)  -  I_temp(n) * I_temp(p) * omega_temp(m) * omega_temp(n)  -  I_temp(m) * I_temp(n) * omega_temp(m) * omega_temp(p) ...
                                    +  I_temp(n) * I_temp(p) * omega_temp(m) * omega_temp(p)  +  I_temp(m) * I_temp(n) * omega_temp(n) * omega_temp(p)  -  I_temp(m) * I_temp(p) * omega_temp(n) * omega_temp(p));
                    [deviazione(test), R_arr(test), L_arr(test)] = motor_parameters(U_temp, omega_temp, I_temp, Ke_arr(test), n0);
                    Ke_mean = Ke_mean + Ke_arr(test);
                    N_test = N_test+1;
                end
            end
        end
        % test with mean of two data set
        Ke_mean = Ke_mean/N_test;
        test = test+1;
        Ke_arr(test) = Ke_mean;
        [deviazione(test), R_arr(test), L_arr(test)] = motor_parameters(U_temp, omega_temp, I_temp, Ke_arr(test), n0);
    end
    
    % check the minimum variance test
    [~, min_dev] = min(deviazione);
    Ke = Ke_arr(min_dev);
    R = R_arr(min_dev);
    L = L_arr(min_dev);
    
    % calc MvK and kL
    if (size(U0,1)>=2) % two set of noload data
        if (omega0(1)~=0)&&(omega0(2)~=0) %both noload data have omega
            MvK = -Ke *( I0(1)*omega0(2) - I0(2)*omega0(1) )/( omega0(1) - omega0(2) );
            kL = -Ke *( I0(2) - I0(1) )/( omega0(1) - omega0(2) );
        else % at least one noload data doesn't have omega
            if (omega0(1)~=0)
                omega0_hat = omega0(1);
            else
                if (omega0(2)~=0)
                    omega0_hat = omega0(2);
                else
                    omega0_hat = 0.985*U0(1)/Ke;
                end
            end
            
            MvK = Ke*( U0(1)*I0(2) - U0(2)*I0(1) )/ ( U0(1) - U0(2) + (R+L*omega0_hat)*(I0(2)-I0(1))  ) ;
            kL = Ke^2*( I0(1) - I0(2) )/ ( U0(1) - U0(2) + (R+L*omega0_hat)*(I0(2)-I0(1))  ) ;
        end
    else
        if (omega0(1)~=0)
            omega0_hat = omega0(1);
        else
            omega0_hat = 0.985*U0(1)/Ke;
        end
        MvK = 0;
        kL = I0(1)*Ke/omega0_hat;
    end
    
    
end

etaG = 1;

% Program
Z = sqrt(R^2+(L*omega)^2);
I = (Mout/etaG + kL*omega + MvK)/Ke;
U = Ke*omega+Z*I;

Pin = U*I;
Pout = Mout*omega;
eta = Pout/Pin;

% check consistency
if (eta<0)||(Pin<0)||(Pout<0)||(U<0)||(I<0)||(Z<0)
    U = NaN;
    I = NaN;
    eta = NaN;
end
    
