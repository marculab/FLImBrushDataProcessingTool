% code to replot interplated video
function [] = replotInterpImage(matName, videoName, radius, alpha, SNR_low)
tic
% read image data
% D= v.Duration-0.04;
D = 0;
v = VideoReader(videoName,'CurrentTime',D);
im = readFrame(v);
%imshow(im)

mat_file = load( matName );
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
assert(v.NumFrames >= max(posData.frames));

% read lifetime related values from mat file
% write into one object
% -lt {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
% -int {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
% -snr {4x1cell}: {[Nx1 double], [Nx1 double], [Nx1 double], [Nx1 double]}
ltData = {};
ltData.lt  = {mat_data(:, 12), mat_data(:, 13), mat_data(:, 14), mat_data(:, 15)};
ltData.it  = {mat_data(:, 16), mat_data(:, 17), mat_data(:, 18), mat_data(:, 19)};
ltData.snr = {mat_data(:, 20), mat_data(:, 21), mat_data(:, 22), mat_data(:, 23)};

[img,scale]=processImg(im, posData, ltData,radius, alpha, SNR_low);

[~,name,~] = fileparts(matName);

for i=1:3
    figure
    imshow(img{i})
    title([name ' Channel', num2str(i)],'Interpreter','none');
    colormap(jet);
    caxis(gca,scale{i});
    h0 = colorbar;
    %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Lifetime (ns)';
    set(gca,'FontSize',20)
    set(gca,'LooseInset',get(gca,'TightInset'))
    saveas(gcf, [name '_ch', num2str(i),'.jpg']);
end
toc
end