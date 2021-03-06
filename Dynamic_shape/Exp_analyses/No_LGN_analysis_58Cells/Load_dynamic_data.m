%%% Import data from spreadsheet
%Workbook: /Users/srivasp/Documents/Work/Crick_2015/Spindle_orientation/Exp_n_model/Dynamic_shape/spindle_dna_cell_data.xlsx
%Worksheet: Sheet1
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.
% Auto-generated by MATLAB on 2018/01/12 16:43:05

%% Import the data
[~, ~, raw] = xlsread('spindle_dna_cell_dynamic_data_jan2018.xlsx','Sheet1');
raw = raw(2:174,4:end);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
Dynamic_shape = table;

%% Allocate imported array to column variable names
Dynamic_shape.id = data(:,1);
Dynamic_shape.frame = data(:,2);
Dynamic_shape.time = data(:,3);
Dynamic_shape.centrosome_a_x = data(:,4);
Dynamic_shape.centrosome_a_y = data(:,5);
Dynamic_shape.centrosome_a_vel = data(:,6);
Dynamic_shape.centrosome_a_dir = data(:,7);
Dynamic_shape.centrosome_b_x = data(:,8);
Dynamic_shape.centrosome_b_y = data(:,9);
Dynamic_shape.centrosome_b_vel = data(:,10);
Dynamic_shape.centrosome_b_dir = data(:,11);
Dynamic_shape.centrosome_mid_x = data(:,12);
Dynamic_shape.centrosome_mid_y = data(:,13);
Dynamic_shape.centrosome_mid_vel = data(:,14);
Dynamic_shape.centrosome_mid_dir = data(:,15);
Dynamic_shape.centrosome_distance = data(:,16);
Dynamic_shape.centrosome_angle = data(:,17);
Dynamic_shape.centrosome_elong = data(:,18);
Dynamic_shape.centrosome_rot = data(:,19);
Dynamic_shape.spindle_pole_a_x = data(:,20);
Dynamic_shape.spindle_pole_a_y = data(:,21);
Dynamic_shape.spindle_pole_a_vel = data(:,22);
Dynamic_shape.spindle_pole_a_dir = data(:,23);
Dynamic_shape.spindle_pole_b_x = data(:,24);
Dynamic_shape.spindle_pole_b_y = data(:,25);
Dynamic_shape.spindle_pole_b_vel = data(:,26);
Dynamic_shape.spindle_pole_b_dir = data(:,27);
Dynamic_shape.spindle_mid_x = data(:,28);
Dynamic_shape.spindle_mid_y = data(:,29);
Dynamic_shape.spindle_mid_vel = data(:,30);
Dynamic_shape.spindle_mid_dir = data(:,31);
Dynamic_shape.spindle_distance = data(:,32);
Dynamic_shape.spindle_angle = data(:,33);
Dynamic_shape.spindle_elong = data(:,34);
Dynamic_shape.spindle_rot = data(:,35);
Dynamic_shape.dna_area = data(:,36);
Dynamic_shape.dna_congression = data(:,37);
Dynamic_shape.dna_x = data(:,38);
Dynamic_shape.dna_y = data(:,39);
Dynamic_shape.dna_length = data(:,40);
Dynamic_shape.dna_width = data(:,41);
Dynamic_shape.dna_angle = data(:,42);
Dynamic_shape.cell_area = data(:,43);
Dynamic_shape.cell_rounding = data(:,44);
Dynamic_shape.cell_x = data(:,45);
Dynamic_shape.cell_y = data(:,46);
Dynamic_shape.cell_length = data(:,47);
Dynamic_shape.cell_shorten = data(:,48);
Dynamic_shape.cell_width = data(:,49);
Dynamic_shape.cell_shrink = data(:,50);
Dynamic_shape.cell_angle = data(:,51);
Dynamic_shape.centrosome_mid_cell_offset = data(:,52);
Dynamic_shape.dna_cell_offset = data(:,53);
Dynamic_shape.centrosome_mid_dna_offset = data(:,54);
Dynamic_shape.centrosome_pole_offset = data(:,55);
Dynamic_shape.cell_spindle_angle_delta = data(:,56);
Dynamic_shape.cell_dna_angle_delta = data(:,57);
Dynamic_shape.spindle_dna_angle_delta = data(:,58);
Dynamic_shape.initial_cell_spindle_angle_delta = data(:,59);
Dynamic_shape.initial_cell_dna_angle_delta = data(:,60);


%% Clear temporary variables
clearvars data raw R;










































