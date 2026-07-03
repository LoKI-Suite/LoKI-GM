% Script to plot two different EEDFs, obtained with LoKI-B, in the same MATLAB graph

% obtain first EEDF
uiwait(msgbox('Select first EEDF file.'));
[file1, path1] = uigetfile('*.txt');
while ~isequal(file1, 'eedf.txt')
  uiwait(warndlg('Please select an EEDF file (''eedf.txt'')','Warning'));
  [file1, path1] = uigetfile('*.txt');
end
fid1 = fopen([path1 file1], 'r');
fgetl(fid1);
data1 = (fscanf(fid1, '%f', [3 inf]));
fclose(fid1);
legend1 = inputdlg('Define legend for first EEDF:');

% obtain second EEDF
uiwait(msgbox('Select second EEDF file.'));
[file2, path2] = uigetfile('*.txt');
while ~isequal(file2, 'eedf.txt')
  uiwait(warndlg('Please select an EEDF file (''eedf.txt'')','Warning'));
  [file2, path2] = uigetfile('*.txt');
end
fid2 = fopen([path2 file2], 'r');
fgetl(fid2);
data2 = (fscanf(fid2, '%f', [3 inf]));
fclose(fid2);
legend2 = inputdlg('Define legend for second EEDF:');

% plot both EEDFs
figure;
semilogy(data1(1,:), data1(2,:), 'r-', data2(1,:), data2(2,:), 'b-');
xlabel('Energy (eV)');
ylabel('Distribution Function (eV^{-3/2})');
legend(gca, {legend1{1} legend2{1}});