[x,y]=meshgrid(linspace(-50,50,1000),linspace(-50,50,1000)); %arbitrary values, just for the simulation
func=zeros(1000,1000); %function for the ring
for n=1:1000
     for m=1:1000
          if sqrt(x(n,m).^2+y(n,m).^2) > 1.8 && sqrt(x(n,m).^2+y(n,m).^2) < 2
              func(n,m)=1;
          end
     end
end 

phase1=exp(1i*1*atan2(y,x)); % turn to 0 if you want only the mode that is created for right circular polarization
phase2=exp(1i*-1*atan2(y,x));  % turn to 0 if you want only the mode created for left circular polarization
func_phase=func.*phase1 + func.*phase2; % this is the field to be scattered out, + is x-polarized, - is y-polarized, circular polarization is if phase1 or phase2 is zero
x_func_phase=func_phase.*cos(atan2(y,x)); % this is the scattered field in the x polarization
y_func_phase=func_phase.*sin(atan2(y,x)); % this is the scattered field in the y polarization

%%
figure(3);
p=sqrt(abs(fftshift(fft2(x_func_phase))).^2+abs(fftshift(fft2(y_func_phase))).^2);
imagesc(sqrt(abs(fftshift(fft2(x_func_phase))).^2));
imagesc(sqrt(abs(fftshift(fft2(y_func_phase))).^2));

axis equal; % the image you should see in fourier
xlim([400 600]);
ylim([400 600]);
colormap("parula")
% colormap("magma")

