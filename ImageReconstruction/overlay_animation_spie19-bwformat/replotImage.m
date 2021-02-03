function [] = replotImage(filename1, filename2)

% read image data
im = imread( filename1 );

% read segmentation position from mat file
% mat file structure:
% Frame,	ID Txt,	Lifetime 1,	Lifetime 2,	Lifetime 3,	Lifetime 4,	Width,	Height,	Z-Axis,	Unknown,	ID Deconv,	Lifetime 1 Deconv,	Lifetime 2 Deconv,	Lifetime 3 Deconv,	Lifetime 4 Deconv,	Intensity 1 Deconv,	Intensity 2 Deconv,	Intensity 3 Deconv,	Intensity 4 Deconv,	SNR 1 Deconv,	SNR 2 Deconv,	SNR 3 Deconv,	SNR 4 Deconv

% write into one object
% -px [Nx1]
% -py [Nx1]
% -frames [integer N]

mat_file = load( filename2 );
fields = fieldnames(mat_file);
mat_data = getfield(mat_file, fields{1,1});
posData = {};
posData.px = mat_data(:, 7);
posData.py = mat_data(:, 8);
posData.frames = size(mat_data, 1);
posData.radius = mat_data(:, 9);

% read lifetime related values from mat file
% write into one object
% -lt {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
% -int {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
% -snr {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
ltData = {};
ltData.lt  = {mat_data(:, 12), mat_data(:, 13), mat_data(:, 14), mat_data(:, 15)};
ltData.it  = {mat_data(:, 16), mat_data(:, 17), mat_data(:, 18), mat_data(:, 19)};
ltData.snr = {mat_data(:, 20), mat_data(:, 21), mat_data(:, 22), mat_data(:, 23)};

% process
processImg(im, posData, ltData);

end