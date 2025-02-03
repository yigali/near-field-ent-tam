%% sum_signal_load
temp=read(Tiff('J2_X1.tif')); %change the name of the file, the file and the code have to be in the same folder
imagesc(temp);
axis equal;
colorbar;
colormap('hot')
xlim([30,130])
ylim([20,120])
