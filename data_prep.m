% load two data files
raw_data1 = importdata('data\velocity-raw-Al1-01-25.mat');
raw_data2 = importdata('data\velocity-raw-Al1-02-09-92.mat');

% remove columns of zeros at the end of the data
velocity_data1 = raw_data1(:, 1:4538);
velocity_data2 = raw_data2(:, 1:4538);

% save data with columns of zeros removed 
save('data\velocity-clean-Al1-01-25', 'velocity_data1')
save('data\velocity-clean-Al1-02-09', 'velocity_data2')