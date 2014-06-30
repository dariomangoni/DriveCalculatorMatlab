clear all
close all
clc

% Setted by optimizator or by user
motor_id = 850;
prefer_benchmark_over_rough_estimate = true; % TRUE if we want to extract data from database instead of using rough estimation*
extract_ALWAYS_parameters_from_benchmark = false; % TRUE if we always want to extract data from database
%* rough estimation is when we have only no load data and some parameters
%without any on-load benchmark

% Benchmark extraction from database
load('database.mat')
mymotor_row = find(motor_table(:,1)==motor_id);
if (isempty(mymotor_row)) || (length(mymotor_row)>1)
    disp(['Motor ID ', num2str(motor_id),' ambiguous']);
    break
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
                    omega0_hat = 0.985*Uo(1)/Ke;
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



%Input
etaG = 1;

U_range = 10;
% omega_range=(1000:100:10000)*pi/30;
omega_range=(1000:100:8000)*pi/30;

I = zeros(length(omega_range),length(U_range));
Mout = I;
Pout = I;
Pin = I;
eta = I;

for omega_ind = 1:length(omega_range)
    for U_ind = 1:length(U_range)
        I(omega_ind,U_ind) = (U_range(U_ind)-Ke*omega_range(omega_ind))./sqrt(R^2+(L*omega_range(omega_ind)).^2);
        Mout(omega_ind,U_ind) = etaG* ( Ke*I(omega_ind,U_ind) - (kL*omega_range(omega_ind)+MvK) );
        Pout(omega_ind,U_ind) = Mout(omega_ind,U_ind).*omega_range(omega_ind);
        Pin(omega_ind,U_ind) = U_range(U_ind) * I(omega_ind,U_ind);
        eta(omega_ind,U_ind) =  Pout(omega_ind,U_ind)./Pin(omega_ind,U_ind);
        
        if eta(omega_ind,U_ind)<0 || eta(omega_ind,U_ind)>1 || ~isreal(eta(omega_ind,U_ind)) || I(omega_ind,U_ind)<0 || Mout(omega_ind,U_ind)<0
            eta(omega_ind,U_ind) = NaN;
            I(omega_ind,U_ind) = NaN;
            Mout(omega_ind,U_ind) = NaN;
            Pin(omega_ind,U_ind) = NaN;
            Pout(omega_ind,U_ind) = NaN;
        end    
        
    end
end

%% plot
if size(U_range,2)>1 && size(omega_range,2)>1
    figure('Name','eta','NumberTitle','off')
    surf(U_range,omega_range/pi*30,eta,'EdgeColor','none');
    zlim([0 1]);
    xlabel('U [V]');
    ylabel('rpm')
    zlabel('eta')
elseif size(omega_range,2)==1
    figure('Name','eta','NumberTitle','off')
    plot(U_range,eta)
    xlabel('U [V]');
    ylabel('eta')
    annotation('textbox','String',['omega: ',num2str(omega_range*30/pi),' rpm'],'Position',[0.5 0.8 0.1 0.1],'FitBoxToText','on');
else
    figure('Name','eta','NumberTitle','off')
    plot(omega_range/pi*30,eta)
    xlabel('rpm');
    ylabel('eta')
    annotation('textbox','String',['U: ',num2str(U_range),' V'],'Position',[0.5 0.8 0.1 0.1],'FitBoxToText','on');
end

%eta max
B = @(w) -Ke*w/(R^2+(L*w)^2);
A = @(w) R^2+(L*w)^2;
C = @(w) (kL*w+MvK);
etamax1 = @(w) ( -(B(w)*Ke+C(w)) + sqrt((B(w)*Ke+C(w)).^2-Ke*B(w).*(B(w)*Ke+C(w))) )/(-Ke/A(w));
etamax2 = @(w) ( -(B(w)*Ke+C(w)) - sqrt((B(w)*Ke+C(w)).^2-Ke*B(w).*(B(w)*Ke+C(w))) )/(-Ke/A(w));



% %%%%%%%%%%%%%%%%%%%%%
% U_range = 5:0.5:50;
% omega_range=(1000:100:10000)'*pi/30;
% 
% I = zeros(length(omega_range),length(U_range));
% Mout = I;
% Pout = I;
% Pin = I;
% eta = I;
% % Program
% for U_ind = 1:length(U_range)
%     I(:,U_ind) = (U_range(U_ind)-Ke*omega_range)./sqrt(R^2+(L*omega_range).^2);
%     Mout(:,U_ind) = etaG* ( Ke*I(:,U_ind) - (kL*omega_range+MvK) );
%     Pout(:,U_ind) = Mout(:,U_ind).*omega_range;
%     Pin(:,U_ind) = U_range(U_ind) * I(:,U_ind);
%     eta(:,U_ind) =  Pout(:,U_ind)./Pin(:,U_ind);
% end
% 
% % Clean up tool
% for m = 1 : size(eta,1)
%     for n = 1 : size(eta,2)
%         if eta(m,n)<0 || eta(m,n)>1 || ~isreal(eta(m,n)) || I(m,n)<0 || Mout(m,n)<0
%             eta(m,n) = NaN;
%             I(m,n) = NaN;
%             Mout(m,n) = NaN;
%             Pin(m,n) = NaN;
%             Pout(m,n) = NaN;
%         end       
%     end
% end
% 
% surf(U_range,omega_range/pi*30,eta,'EdgeColor','none');
% zlim([0 1]);
% xlabel('U [V]');
% ylabel('rpm')
% zlabel('eta')
