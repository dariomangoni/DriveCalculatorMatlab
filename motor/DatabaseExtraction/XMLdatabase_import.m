close all
clear all
clc

% XML DOM "doubles" the entries of the database, but I don't know why.
% It's not a problem due to CarryReturn+LineFeed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Motors parameters acquisition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmldoc = xmlread('XMLdatabase\MotorsDatabase.xml');
element_list = getDocumentElement(xmldoc); % element_list now has the motor list
element_total = getLength(element_list); % element_total is the total number of motors in the database

motor_table = zeros(fix(element_total/2),19);
row_index = 0;
for element_index = 1:element_total
    if hasChildNodes(element_list.item(element_index-1))
        row_index = row_index+1;
        element_data_list = getChildNodes(item( element_list, element_index-1 )); %N.B. element_list has 0 as first index; element_data_list contains my_id, tbs, Name, ecc
        element_data_total = getLength(element_data_list); %number of data of each motor (doubled)
        
        for element_data_index = 1:element_data_total %we have only 19 real entries (from myid to ts, excluding Name and Rem)
            if hasChildNodes(element_data_list.item(element_data_index-1))
                element_data = item( element_data_list, element_data_index-1 );
                data_name = char(getNodeName(element_data));
                switch data_name
                    case 'myid'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,1) = str2double(data_text);
                    case 'tbs'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,2) = str2double(data_text);
                    case 'custom'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,3) = str2double(data_text);
                    case 'meas'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,4) = str2double(data_text);
                    case 'match'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,5) = str2double(data_text);
                    case 'StatDia'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,6) = str2double(data_text);
                    case 'StatH'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,7) = str2double(data_text);
                    case 'Turns'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,8) = str2double(data_text);
                    case 'WireDia'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,9) = str2double(data_text);
                    case 'Delta'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,10) = str2double(data_text);
                    case 'Imax'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,11) = str2double(data_text);
                    case 'Weight'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,12) = str2double(data_text);
                    case 'mgear_id'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,13) = str2double(data_text);
                    case 'mesc_id'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,14) = str2double(data_text);
                    case 'kV'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,15) = str2double(data_text);
                    case 'Rm'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,16) = str2double(data_text);
                    case 'K'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,17) = str2double(data_text);
                    case 'MvK'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,18) = str2double(data_text);
                    case 'kL'
                        data_text = getTextContent(element_data);
                        motor_table(row_index,19) = str2double(data_text);
                end
            end

        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Motor benchmarks acquisition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmldoc = xmlread('XMLdatabase\MDataDatabase.xml');
element_list = getDocumentElement(xmldoc); % element_lsit now has the motor list
element_total = getLength(element_list); % element_total is the total number of motors in the database

mdata_table = zeros(fix(element_total/2),10);
row_index = 0;
for element_index = 1:element_total
    if hasChildNodes(element_list.item(element_index-1))
        row_index = row_index+1;
        element_data_list = getChildNodes(item( element_list, element_index-1 )); %N.B. element_list has 0 as first index; element_data_list contains my_id, tbs, Name, ecc
        element_data_total = getLength(element_data_list); %number of data of each mdata (it must be 19)
        for element_data_index = 1:element_data_total
            if hasChildNodes(element_data_list.item(element_data_index-1))
                element_data = item( element_data_list, element_data_index-1 );
                data_name = char(getNodeName(element_data));
                switch data_name
                    case 'myid'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,1) = str2double(data_text);
                    case 'tbs'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,2) = str2double(data_text);
                    case 'motor_id'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,3) = str2double(data_text);
                    case 'U'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,4) = str2double(data_text);
                    case 'n'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,5) = str2double(data_text);
                    case 'I'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,6) = str2double(data_text);
                    case 'IsNoLoad'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,7) = str2double(data_text);
                    case 'prop_id'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,8) = str2double(data_text);
                    case 'Alt'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,9) = str2double(data_text);
                    case 'Temp'
                        data_text = getTextContent(element_data);
                        mdata_table(row_index,10) = str2double(data_text);
                end
            end

        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Batteries benchmarks acquisition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmldoc = xmlread('XMLdatabase\BatteriesDatabase.xml');
element_list = getDocumentElement(xmldoc); % element_lsit now has the motor list
element_total = getLength(element_list); % element_total is the total number of motors in the database

batteries_table = zeros(fix(element_total/2),8);
row_index = 0;
for element_index = 1:element_total
    if hasChildNodes(element_list.item(element_index-1))
        row_index = row_index+1;
        element_data_list = getChildNodes(item( element_list, element_index-1 )); %N.B. element_list has 0 as first index; element_data_list contains my_id, tbs, Name, ecc
        element_data_total = getLength(element_data_list); %number of data of each mdata (it must be 19)
        for element_data_index = 1:element_data_total
            if hasChildNodes(element_data_list.item(element_data_index-1))
                element_data = item( element_data_list, element_data_index-1 );
                data_name = char(getNodeName(element_data));
                switch data_name
                    case 'myid'
                        data_text = getTextContent(element_data);
                        batteries_table(row_index,1) = str2double(data_text);
                    case 'tbs'
                        data_text = getTextContent(element_data);
                        batteries_table(row_index,2) = str2double(data_text);
                    case 'Imax'
                        data_text = getTextContent(element_data);
                        batteries_table(row_index,3) = str2double(data_text);
                    case 'Ipeak'
                        data_text = getTextContent(element_data);
                        batteries_table(row_index,4) = str2double(data_text);
                    case 'Capacity'
                        data_text = getTextContent(element_data);
                        batteries_table(row_index,5) = str2double(data_text);
                    case 'Weight'
                        data_text = getTextContent(element_data);
                        batteries_table(row_index,6) = str2double(data_text);
                    case 'Ri'
                        data_text = getTextContent(element_data);
                        batteries_table(row_index,7) = str2double(data_text);
                    case 'Volt'
                        data_text = getTextContent(element_data);
                        batteries_table(row_index,8) = str2double(data_text);
                end
            end

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DriveSets benchmarks acquisition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmldoc = xmlread('XMLdatabase\DriveSetsDatabase.xml');
element_list = getDocumentElement(xmldoc); % element_lsit now has the motor list
element_total = getLength(element_list); % element_total is the total number of motors in the database

drivesets_table = zeros(fix(element_total/2),8);
row_index = 0;
for element_index = 1:element_total
    if hasChildNodes(element_list.item(element_index-1))
        row_index = row_index+1;
        element_data_list = getChildNodes(item( element_list, element_index-1 )); %N.B. element_list has 0 as first index; element_data_list contains my_id, tbs, Name, ecc
        element_data_total = getLength(element_data_list); %number of data of each mdata (it must be 19)
        for element_data_index = 1:element_data_total
            if hasChildNodes(element_data_list.item(element_data_index-1))
                element_data = item( element_data_list, element_data_index-1 );
                data_name = char(getNodeName(element_data));
                switch data_name
                    case 'myid'
                        data_text = getTextContent(element_data);
                        drivesets_table(row_index,1) = str2double(data_text);
                    case 'tbs'
                        data_text = getTextContent(element_data);
                        drivesets_table(row_index,2) = str2double(data_text);
                    case 'esc_id'
                        data_text = getTextContent(element_data);
                        drivesets_table(row_index,3) = str2double(data_text);
                    case 'bat_id'
                        data_text = getTextContent(element_data);
                        drivesets_table(row_index,4) = str2double(data_text);
                    case 'mot_id'
                        data_text = getTextContent(element_data);
                        drivesets_table(row_index,5) = str2double(data_text);
                    case 'gear_id'
                        data_text = getTextContent(element_data);
                        drivesets_table(row_index,6) = str2double(data_text);
                    case 'prop_id'
                        data_text = getTextContent(element_data);
                        drivesets_table(row_index,7) = str2double(data_text);
                    case 'last_U'
                        data_text = getTextContent(element_data);
                        drivesets_table(row_index,8) = str2double(data_text);
                end
            end

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESC benchmarks acquisition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmldoc = xmlread('XMLdatabase\ESCDatabase.xml');
element_list = getDocumentElement(xmldoc); % element_lsit now has the motor list
element_total = getLength(element_list); % element_total is the total number of motors in the database

ESC_table = zeros(fix(element_total/2),8);
row_index = 0;
for element_index = 1:element_total
    if hasChildNodes(element_list.item(element_index-1))
        row_index = row_index+1;
        element_data_list = getChildNodes(item( element_list, element_index-1 )); %N.B. element_list has 0 as first index; element_data_list contains my_id, tbs, Name, ecc
        element_data_total = getLength(element_data_list); %number of data of each mdata (it must be 19)
        for element_data_index = 1:element_data_total
            if hasChildNodes(element_data_list.item(element_data_index-1))
                element_data = item( element_data_list, element_data_index-1 );
                data_name = char(getNodeName(element_data));
                switch data_name
                    case 'myid'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,1) = str2double(data_text);
                    case 'tbs'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,2) = str2double(data_text);
                    case 'Imax'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,3) = str2double(data_text);
                    case 'Ipeak'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,4) = str2double(data_text);
                    case 'R'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,5) = str2double(data_text);
                    case 'Weight'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,6) = str2double(data_text);
                    case 'NiMHmin'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,7) = str2double(data_text);
                    case 'NiMHmax'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,8) = str2double(data_text);
                    case 'LiPomin'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,9) = str2double(data_text);
                    case 'LiPomax'
                        data_text = getTextContent(element_data);
                        ESC_table(row_index,10) = str2double(data_text);
                end
            end

        end
    end
end

%% Salvataggio dati
save('database.mat','mdata_table', 'motor_table', 'batteries_table', 'drivesets_table', 'ESC_table');