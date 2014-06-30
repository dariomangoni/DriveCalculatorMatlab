clear all
clc

global motor_table %#ok<*NUSED>
load('database.mat')

motor_total = size(motor_table,1); % total number of motors

omega = 3000*pi/30;
Mout = 0.15;

motor_index = 0;
for motor_id = 1:motor_total
    motor_index = motor_index+1;
    [U(motor_index), I(motor_index), eta(motor_index)] = motor_function( omega, Mout, motor_id );
    motor_id_list(motor_index) = motor_id;
    if isnan(U(motor_index))
        motor_index = motor_index-1;
    end
end

%% plot
close all
figure('Name','eta','NumberTitle','off')
plot(motor_id_list,eta,'gx')
figure('Name','U','NumberTitle','off')
plot(motor_id_list,U,'bx')
figure('Name','I','NumberTitle','off')
plot(motor_id_list,I,'rx')
figure('Name','Power','NumberTitle','off')
plot(motor_id_list,U.*I,'kx')

comparison_table = sortrows([motor_id_list', U', I', eta'],4);
