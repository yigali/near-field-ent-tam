warning('off','all')
warning
%% properties
image_num = 86400;
x_size = 46;
y_size = 46;
% x_sizeS = 46;
% y_sizeS = 61;
% SignalA=zeros(x_sizeS,y_sizeS,image_num);
Signal=zeros(x_size,y_size,image_num);

cor_matS_A_batch = zeros(x_size, y_size, x_size, y_size);
cor_matS_B_batch = zeros(x_size, y_size, x_size, y_size);
cor_matS_A = zeros(x_size, y_size, x_size, y_size);
cor_matS_B = zeros(x_size, y_size, x_size, y_size);

% %% calculate corrolations
Batch = {'A_X','B_X','C_X'};
% yRange = 2:42; % Rows 2 to 42 (inclusive)
% xRange = 8:48; % Columns 8 to 48 (inclusive)

for r=1:length(Batch)
    for l=1:image_num
        disp(l)
        temp=read(Tiff([char(Batch(1,r)), num2str(l) '.tif']));
        Signal(:,:,l)=double(temp);
%         Signal(:,:,l)=SignalA(yRange, xRange,l);
    end
    for i=1:x_size
        count = ['i =' ,num2str(i)];
        disp(count); 
        for j=1:y_size
            alphaS=zeros(x_size, y_size);
            betaS=zeros(x_size, y_size);
                parfor k=1:(image_num-1)
                    alphaS=alphaS+Signal(i,j,k)*Signal(:,:,k);
                    betaS=betaS+Signal(i,j,k)*Signal(:,:,k+1); 
                end
            cor_matS_A_batch(i,j,:,:) = alphaS/(image_num-1);
            cor_matS_B_batch(i,j,:,:) = betaS/(image_num-2);

            cor_matS_A(i,j,:,:)=cor_matS_A(i,j,:,:)+cor_matS_A_batch(i,j,:,:);
            cor_matS_B(i,j,:,:)=cor_matS_B(i,j,:,:)+cor_matS_B_batch(i,j,:,:);
        end
    end
end

% calculate the final JPD corrolations
cor_matS_C = zeros(x_size, y_size, x_size, y_size);
for i=1:x_size
    disp(i)
    for j=1:y_size
        cor_matS_C(i,j,:,:)=cor_matS_A(i,j,:,:)-cor_matS_B(i,j,:,:);
    end
end