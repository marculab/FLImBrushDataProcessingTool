% code to replot interplated video
function [] = replotInterpImage(filename1, filename2)
tic
% read image data
v = VideoReader(filename1);
D= v.Duration-0.04;
v = VideoReader(filename1,'CurrentTime',D);
im = readFrame(v);
imshow(im)

mat_file = load( filename2 );
fields = fieldnames(mat_file);
mat_data = getfield(mat_file, fields{1,1});
posData = {};
posData.px = mat_data(:, 7);
posData.py = mat_data(:, 8);
posData.frames = mat_data(:,1);
posData.radius = mat_data(:,9);
% align dimensions of segmentation positions with video
% It seems that the txtfile has same dimension with video frame
% number; awaiting confirmation

%assert(v.NumberOfFrames == posData.frames);
assert(v.NumberOfFrames >= max(posData.frames));

% read lifetime related values from mat file
% write into one object
% -lt {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
% -int {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
% -snr {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
ltData = {};
ltData.lt  = {mat_data(:, 12), mat_data(:, 13), mat_data(:, 14), mat_data(:, 15)};
ltData.it  = {mat_data(:, 16), mat_data(:, 17), mat_data(:, 18), mat_data(:, 19)};
ltData.snr = {mat_data(:, 20), mat_data(:, 21), mat_data(:, 22), mat_data(:, 23)};

processImg(im, posData, ltData, filename2);

toc
end