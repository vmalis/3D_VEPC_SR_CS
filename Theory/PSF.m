% PSF for checking incoherence


N = 64;                    % Size of a sunction
D = N/2;                    % to indicate origin at the center of the function
w = 3;                    % decay rate for exponential function
y = repmat(1:N,N,1);
x = y';
r = sqrt((20-x).^2+(5-y).^2);    % definition of radius

psf = (airy(w*(r)/2).^2)/pi;
psf = psf/max(psf(:));

figure;
set(gcf,'color','w');
mesh(x,y,psf,'EdgeColor','black')
%surf(x,y,psf,'FaceColor','interp','FaceLighting','gouraud')
axis tight
axis square
set(gca,'xtick',0:8:64)
ztickformat('%,.2f')
set(gca,'ytick',0:8:64)
set(gca,'ztick',-1:0.25:1)
zlim([-1,1])
xlim([0,64])
ylim([0,64])
axis off
export_fig('1.eps')


XFM = Wavelet('Daubechies',4,4);

W=XFM'*psf;
W=W./max(W(:));
figure;
set(gcf,'color','w');
mesh(x,y,W,'EdgeColor','black')
%surf(x,y,W,'FaceColor','interp','FaceLighting','gouraud')
axis tight
axis square
set(gca,'xtick',0:8:64)
ztickformat('%,.2f')
set(gca,'ytick',0:8:64)
set(gca,'ztick',-1:0.25:1)
zlim([-1,1])
xlim([0,64])
ylim([0,64])
axis off
export_fig('2.eps')



Q=fftshift(fftshift(fft2(W),1),2);
figure;
set(gcf,'color','w');
Q=Q./max(Q(:));
mesh(x,y,abs(Q),'EdgeColor','black')
%surf(x,y,real(Q),'FaceColor','interp','FaceLighting','gouraud')
axis tight
axis square
set(gca,'xtick',0:8:64)
ztickformat('%,.2f')
set(gca,'ytick',0:8:64)
set(gca,'ztick',-1:0.25:1)
zlim([-1,1])
xlim([0,64])
ylim([0,64])
axis off
export_fig('3.eps')




M = randi([0 1], 64,64);
temp=M+1;
temp(temp==2)=0;
C=repmat(temp,[1,1,2]);
C=cat(3,ones(size(M)),C);


figure;
set(gcf,'color','w');
Q=Q./max(Q(:));
surf(x,y,abs(Q),C,'EdgeColor','black')
%surf(x,y,real(Q),'FaceColor','interp','FaceLighting','gouraud')
axis tight
axis square
set(gca,'xtick',0:8:64)
ztickformat('%,.2f')
set(gca,'ytick',0:8:64)
set(gca,'ztick',-1:0.25:1)
zlim([-1,1])
xlim([0,64])
ylim([0,64])
axis off
export_fig('4.eps')





Q_us=Q.*M;
W_us=ifft2c(Q_us);
W_us=W_us./max(W_us(:));

figure;
set(gcf,'color','w');
mesh(x,y,real(flip(flip(W_us',1),2)),'EdgeColor','black')
axis tight
axis square
set(gca,'xtick',0:8:64)
ztickformat('%,.2f')
set(gca,'ytick',0:8:64)
set(gca,'ztick',-1:0.25:1)
zlim([-1,1])
xlim([0,64])
ylim([0,64])
axis off
export_fig('5.eps')


psf_us=XFM*W_us;
psf_us=psf_us./max(psf_us(:));
figure;
set(gcf,'color','w');
mesh(x,y,real(psf_us),'EdgeColor','black')
axis tight
axis square
set(gca,'xtick',0:8:64)
ztickformat('%,.2f')
set(gca,'ytick',0:8:64)
set(gca,'ztick',-1:0.25:1)
zlim([-1,1])
xlim([0,64])
ylim([0,64])
axis off
export_fig('6.eps')



