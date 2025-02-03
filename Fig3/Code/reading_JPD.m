warning('off')
clc;
clear all;

%% filter
%load the matrixs befor ranning this code
%crate filter for the plasmons scattering
fil=zeros(26,26);
for l=1:26
    for m=1:26
        if (l-12.5)^2 +(m-11)^2 > 6^2 && (l-12.5)^2 +(m-11)^2 < 9.5^2
            fil(l,m)=1;
        end
    end
end
figure();
imagesc(fil);
hold on
axis equal
xlim([0,22])
viscircles([11 12.5],6,"LineStyle","--","Color","r","LineWidth",0.1);
viscircles([11 12.5],10,"LineStyle","--","Color","r","LineWidth",0.1);
ylim([2,24])

%%        
% relevent pixels from JPD for the video
pixel_mat_x = [17, 18, 19, 20, 20, 20, 20, 20, 20, 20, 21, 20, 19, 18, 17, 17, 16, 15, 14, 13, 12, 10, 9, 8, 7, 6, 4, 3, 4, 5, 6, 6, 6, 6, 6, 7, 8, 9, 10, 11, 11, 12, 13, 13, 14, 13,17];
pixel_mat_y = [19, 17, 16, 15, 14, 13, 12, 11, 10, 9, 10, 8, 6, 5, 4, 5, 4, 4, 3, 3, 3, 3, 3, 4, 4, 6, 9, 9, 10, 11, 12, 14, 15, 16, 17, 15, 17, 18, 18, 20, 18, 18, 18, 19, 20, 20,19];
num_of_pixels = length(pixel_mat_y);

sum_mat = zeros(26, 26);

for n=1:num_of_pixels
    
    x_coord = pixel_mat_x(n);
    y_coord = pixel_mat_y(n);
    % chose linear or circuler matrix from JPD
    y_original=reshape(cor_matS_C_linear(x_coord,y_coord,:,:), 26, 26); % linear 
%     y_original=reshape(cor_matS_C_circular(x_coord,y_coord,:,:), 26, 26); %circular
    y = y_original;
    
    % zero out pixels near the origin
    for i=-2:1:5
        for j=-2:1:5
            y(x_coord+i, y_coord+j) = 1;

        end
    end

    z=real(log(y));
    t=sign(y);
            
    % zero out mask with negative entries
    for i=1:26

        for j=1:26
                    
            if sign(t(i, j)) < 0
                t(i, j) = 0;
            end
        end
    end


    sum_mat = sum_mat + z.*t.*fil;
    file_name = sprintf('%d.tif', n);
    full_file_path = fullfile('C:\Users\amkam\OneDrive\שולחן העבודה\מדידות\linear', file_name); 
    
    imagesc(sum_mat)
    axis equal
    xlim([0,22])
    ylim([2,24])
    colorbar()
    colormap('hot')
    caxis([0,5])
    hold all
    plot(pixel_mat_y(n),pixel_mat_x(n), '-p','MarkerFaceColor','green', 'MarkerSize',15);
    viscircles([11 12.5],6,"LineStyle","--","Color","r","LineWidth",0.1);
    viscircles([11 12.5],10,"LineStyle","--","Color","r","LineWidth",0.1);
    % Write the matrix to a TIFF image
    saveas(gcf, full_file_path, 'tif');

end


% Make sure pictures names are only numbers, and that their numbered at the right order!
video = VideoWriter('video_name1.mp4','MPEG-4');  %format
num_pics=num_of_pixels; %number of pictures
video.FrameRate=4; %variable - play with it until it looks good
open(video);
for s=1:num_pics
    file_name = strcat('C:\Users\amkam\OneDrive\שולחן העבודה\מדידות\linear\', sprintf('%d.tif', s)); %change the location to save
    I = imread(file_name);
    writeVideo(video,I); 
end
close(video); 





