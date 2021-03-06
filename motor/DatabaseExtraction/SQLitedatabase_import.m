close all
clear all
clc

%% SQLite database import
addpath('mksqlite-1.13')
mksqlite('open', 'SQLitedatabase\DCbase.dcd');
mksqlite('PRAGMA synchronous = 2'); % it ensures that also with a computer crash database won't be corrupted

motor_struct     = mksqlite('SELECT * FROM Motors');
batteries_struct = mksqlite('SELECT * FROM Batteries');
mdata_struct     = mksqlite('SELECT * FROM MData');
ESC_struct       = mksqlite('SELECT * FROM ESC');
drivesets_struct = mksqlite('SELECT * FROM DriveSets');

mksqlite('close');

%% Extract tables from structures

motor_table = zeros(size(motor_struct,1),19);
for index_count = 1:size(motor_struct,1)
    motor_table(index_count,:) = [ motor_struct(index_count).myid,...
                                   motor_struct(index_count).tbs,...
                                   motor_struct(index_count).custom,...
                                   motor_struct(index_count).meas,...
                                   motor_struct(index_count).match,...
                                   motor_struct(index_count).StatDia,...
                                   motor_struct(index_count).StatH,...
                                   motor_struct(index_count).Turns,...
                                   motor_struct(index_count).WireDia,...
                                   motor_struct(index_count).Delta,...
                                   motor_struct(index_count).Imax,...
                                   motor_struct(index_count).Weight,...
                                   motor_struct(index_count).mgear_id,...
                                   motor_struct(index_count).mesc_id,...
                                   motor_struct(index_count).kV,...
                                   motor_struct(index_count).Rm,...
                                   motor_struct(index_count).K,...
                                   motor_struct(index_count).MvK,...
                                   motor_struct(index_count).kL ];
end

mdata_table = zeros(size(mdata_struct,1),10);
for index_count = 1:size(mdata_struct,1)
    mdata_table(index_count,:) = [ mdata_struct(index_count).myid,...
                                   mdata_struct(index_count).tbs,...
                                   mdata_struct(index_count).motor_id,...
                                   mdata_struct(index_count).U,...
                                   mdata_struct(index_count).n,...
                                   mdata_struct(index_count).I,...
                                   mdata_struct(index_count).IsNoLoad,...
                                   mdata_struct(index_count).prop_id,...
                                   mdata_struct(index_count).Alt,...
                                   mdata_struct(index_count).Temp ];
end

ESC_table = zeros(size(ESC_struct,1),10);
for index_count = 1:size(ESC_struct,1)
    ESC_table(index_count,:) = [ ESC_struct(index_count).myid,...
                                 ESC_struct(index_count).tbs,...
                                 ESC_struct(index_count).Imax,...
                                 ESC_struct(index_count).Ipeak,...
                                 ESC_struct(index_count).R,...
                                 ESC_struct(index_count).Weight,...
                                 ESC_struct(index_count).NiMHmin,...
                                 ESC_struct(index_count).NiMHmax,...
                                 ESC_struct(index_count).LiPomin,...
                                 ESC_struct(index_count).LiPomax ];
end

drivesets_table = zeros(size(drivesets_struct,1),8);
for index_count = 1:size(drivesets_struct,1)
    drivesets_table(index_count,:) = [ drivesets_struct(index_count).myid,...
                                       drivesets_struct(index_count).tbs,...
                                       drivesets_struct(index_count).esc_id,...
                                       drivesets_struct(index_count).bat_id,...
                                       drivesets_struct(index_count).mot_id,...
                                       drivesets_struct(index_count).gear_id,...
                                       drivesets_struct(index_count).prop_id,...
                                       drivesets_struct(index_count).last_U ];
end

batteries_table = zeros(size(batteries_struct,1),8);
for index_count = 1:size(batteries_struct,1)
    batteries_table(index_count,:) = [ batteries_struct(index_count).myid,...
                                       batteries_struct(index_count).tbs,...
                                       batteries_struct(index_count).Imax,...
                                       batteries_struct(index_count).Ipeak,...
                                       batteries_struct(index_count).Capacity,...
                                       batteries_struct(index_count).Weight,...
                                       batteries_struct(index_count).Ri,...
                                       batteries_struct(index_count).Volt ];
end



%% Salvataggio dati
save('database.mat','mdata_table', 'motor_table', 'batteries_table', 'drivesets_table', 'ESC_table');
save('database_struct.mat','mdata_struct', 'motor_struct', 'batteries_struct', 'drivesets_struct', 'ESC_struct');