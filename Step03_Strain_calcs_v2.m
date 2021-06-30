%% -----------------------------------------------------------------------   
%                       STRAIN FOR CS VEPC data
%
%  UC San Diego / March 2019 / Vadim Malis
%
%  description: Top level script 3
%
% ------------------------------------------------------------------------
clc
clear all

%% volume merge
load('force_data.mat')
load('cs_processed.mat')
cd DTI
load('DTI.mat')
cd ..

% important parameters
MVC=unique({Recon.MVC},'stable');
num_iter=20; kappa=20; option=1; delta_t=3/44;
pixelspacing = [1.1719,1.1719,5]';
%num_iter=10; kappa=2; option=1; delta_t=1/7;


for volume=1:size(unique({Recon.MVC}),2)

    Data(volume).ID = unique([Recon.ID],'stable');
    Data(volume).MVC= MVC(volume);
    
    same_mvc_series = find(contains({Recon.MVC},MVC(volume)));
    
    %[~,~,dti_locations] = intersect(round([Recon(same_mvc_series).location]),round(DTI_data.location'),'stable');
    slicDTI = 3;
    dti_locations=[];
    
    for tt=1:slicDTI
        x=abs(repmat(round(Recon(tt).location),[1,size(DTI_data.location,2)])-round(DTI_data.location));
        [~,temp]=min(x);
        dti_locations=[dti_locations,temp];
    end
    
    %dti_locations = [6,7,8];
    
    
    
    Data(volume).vx = permute(cat(4,Recon(same_mvc_series(1)).Vx, Recon(same_mvc_series(1)).Vx,Recon(same_mvc_series(2)).Vx,Recon(same_mvc_series(3)).Vx,Recon(same_mvc_series(3)).Vx),[1,2,4,3]);
    Data(volume).vy = permute(cat(4,Recon(same_mvc_series(1)).Vy, Recon(same_mvc_series(1)).Vy,Recon(same_mvc_series(2)).Vy,Recon(same_mvc_series(3)).Vy,Recon(same_mvc_series(3)).Vy),[1,2,4,3]);
    Data(volume).vz = permute(cat(4,Recon(same_mvc_series(1)).Vz, Recon(same_mvc_series(1)).Vz,Recon(same_mvc_series(2)).Vz,Recon(same_mvc_series(3)).Vz,Recon(same_mvc_series(3)).Vz),[1,2,4,3]);
    
    
    
    
    mask=zeros(size(Data(volume).vx));
    mask(Data(volume).vx~=0)=1;
    Data(volume).mask = mask;
    Data(volume).m  = mask.*permute(cat(4,Recon(same_mvc_series(1)).M, Recon(same_mvc_series(1)).M,Recon(same_mvc_series(2)).M,Recon(same_mvc_series(3)).M,Recon(same_mvc_series(3)).M),[1,2,4,3]);
    
    Data(volume).dti_vec  = DTI_data.DTI_eigenvector(:,:,dti_locations,:,:);
    Data(volume).dti_val  = DTI_data.DTI_eigenvalue(:,:,dti_locations,:);
    
    [~,~,force_series] = intersect({Recon(same_mvc_series).series_num},{Force.series_num},'stable');
    
    Data(volume).MVC_true = mean([Force(force_series).pcent]);
    Data(volume).force    = mean(reshape([Force(force_series).mean],[size(Force(force_series(1)).mean,2),size(force_series,1)]),2);
    
    frames=size(Data(volume).vx,4);
    
    
    h=waitbar(0,'Filtering');
    
    for k=1:frames

        
         Data(volume).vx_sm(:,:,:,k) = anisodiff3D(Data(volume).vx(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);
         Data(volume).vy_sm(:,:,:,k) = anisodiff3D(Data(volume).vy(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);
         Data(volume).vz_sm(:,:,:,k) = anisodiff3D(Data(volume).vz(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);

%         num_iter=10; kappa=2; option=1; delta_t=1/7;
% 
%         for slice=1:3
%             Data(volume).vx_sm(:,:,slice,k) = anisodiff2D(Data(volume).vx(:,:,slice,k),...
%                               num_iter,delta_t,kappa,option);
%             Data(volume).vy_sm(:,:,slice,k) = anisodiff2D(Data(volume).vy(:,:,slice,k),...
%                               num_iter,delta_t,kappa,option);
%             Data(volume).vz_sm(:,:,slice,k) = anisodiff2D(Data(volume).vz(:,:,slice,k),...
%                               num_iter,delta_t,kappa,option);
%         end




         waitbar(k/frames,h)
         
    end
    
    Data(1).vx = Recon(j).Vx_SM
    
    Data(volume).vx(:,:,5,:)=[];
    Data(volume).vx(:,:,1,:)=[];
    Data(volume).vx_sm(:,:,5,:)=[];
    Data(volume).vx_sm(:,:,1,:)=[];
    
    Data(volume).vy(:,:,5,:)=[];
    Data(volume).vy(:,:,1,:)=[];
    Data(volume).vy_sm(:,:,5,:)=[];
    Data(volume).vy_sm(:,:,1,:)=[];
    
    Data(volume).vz(:,:,5,:)=[];
    Data(volume).vz(:,:,1,:)=[];
    Data(volume).vz_sm(:,:,5,:)=[];
    Data(volume).vz_sm(:,:,1,:)=[];
    
    Data(volume).m(:,:,5,:)=[];
    Data(volume).m(:,:,1,:)=[];
    
    Data(volume).mask(:,:,5,:)=[];
    Data(volume).mask(:,:,1,:)=[];
    
    
    close(h)
end

clearvars -except Data frames pixelspacing
close all








%% ROI and calcs


a=figure;
imshow(mat2gray(Data(end).m(:,:,2,1)),'Initialmagnification',400)
h_all = drawfreehand(gca);
mask_all = roiWait(h_all);
close(a)

a=figure;
imshow(mat2gray(Data(end).m(:,:,2,1)),'Initialmagnification',400);
h_gm  = drawrectangle(gca,'Position',[50 50 10 30]);
mask_gm = roiWait(h_gm);
close(a)

a=figure;
imshow(mat2gray(Data(end).m(:,:,2,1)),'Initialmagnification',400);h_all
h_sol  = drawrectangle(gca,'Position',[50 50 10 30]);
mask_sol = roiWait(h_sol);
close(a)


for series=1:size(Data,2)

%predifined from acqusition
trigger=ones(frames-1,1)*0.136;


% create coordinates

mask=Data(series).mask(:,:,:,:);
mask(mask==0)=NaN;
[x,y,z]=meshgrid(1:256,1:256,1:3);

x(isnan(mask(:,:,:,1)))=NaN;
y(isnan(mask(:,:,:,1)))=NaN;
z(isnan(mask(:,:,:,1)))=NaN;


% track each pixel in the loop
% pixel is skipped if nan
% input coordinates are from the 2nd slice out of 3 (middle slice)

XS=zeros([size(x),frames]);
YS=XS;
ZS=XS;
VX=XS;
VY=XS;
VZ=XS;
VR=XS;

    for yi=1:256
        yi
          for xi=1:256
              
              for zi=1:3
              
                [xs,ys,zs,vx,vy,vz,vr]=track3dv2(x(yi,xi,zi),y(yi,xi,zi),z(yi,xi,zi),Data(series).vx_sm,Data(series).vz_sm,Data(series).vy_sm,trigger,pixelspacing);
          
                XS(yi,xi,zi,:)=xs;
                YS(yi,xi,zi,:)=ys;
                ZS(yi,xi,zi,:)=zs;
                VX(yi,xi,zi,:)=vx;
                VY(yi,xi,zi,:)=vy;
                VZ(yi,xi,zi,:)=vz;
                VR(yi,xi,zi,:)=vr;
                
              end
          
          end
    end


mask_sol_dynamic=repmat(mask_sol,[1,1,3,frames]);
mask_gm_dynamic=repmat(mask_gm,[1,1,3,frames]);


%%                             
% tensor calcs                                
Data(series).strain     =   strain_tensor_calcs(VX,Data(series).vx_sm, VY, Data(series).vz_sm, VZ, Data(series).vy_sm,XS,YS,ZS,pixelspacing,trigger);
[Data(series).strain.ED, Data(series).strain.LD, Data(series).strain.SRD]  =   fa_strain_calcs(Data(series).strain,Data(series).dti_vec);

% masking
Data(series).strain_SOL = strain_ROI(Data(series).strain,double(mask_sol_dynamic));
Data(series).strain_GM  = strain_ROI(Data(series).strain,double(mask_gm_dynamic));

% peak
Data(series).strain_SOL_peak = peak_ROI(Data(series).strain_SOL);
Data(series).strain_GM_peak  = peak_ROI(Data(series).strain_GM);


%colormaps
foldername=Data(series).MVC{1};
mkdir(foldername)
cd(foldername)
m=repmat(Data(series).m(:,:,2,1),[1,1,1,size(Data(series).m,4)]);

% displacements

limits=[0,7];
suffix=[foldername,'_dx'];
units = '$\Delta_{x} \quad  [\mathrm{mm}]$';
d=10*squeeze(Data(series).strain.dx(:,:,2,:));
colormap_strain(m,abs(d),suffix,limits,units);

suffix=[foldername,'_dy'];
units = '$\Delta_{y} \quad  [\mathrm{mm}]$';
d=10*squeeze(Data(series).strain.dy(:,:,2,:));
colormap_strain(m,abs(d),suffix,limits,units);

suffix=[foldername,'_dz'];
units = '$\Delta_{z} \quad  [\mathrm{mm}]$';
d=10*squeeze(Data(series).strain.dz(:,:,2,:));
colormap_strain(m,abs(d),suffix,limits,units);


% velocities
limits=[-3,3];
suffix=[foldername,'_vx'];
units = '$v_{x} \quad  [\mathrm{cm/s}]$';
d=squeeze(Data(series).strain.vx(:,:,2,:));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_vy'];
units = '$v_{y} \quad  [\mathrm{cm/s}]$';
d=squeeze(Data(series).strain.vy(:,:,2,:));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_vz'];
units = '$v_{z} \quad  [\mathrm{cm/s}]$';
d=squeeze(Data(series).strain.vz(:,:,2,:));
colormap_strain(m,d,suffix,limits,units);


% Strain Euler
limits=[-0.7,0.7];
suffix=[foldername,'_Eneg'];
units = '$E_{\lambda_{1}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.E_lambda(:,:,:,:,1));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Esum'];
units = '$E_{\lambda_{2}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.E_lambda(:,:,:,:,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Epos'];
units = '$E_{\lambda_{3}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.E_lambda(:,:,:,:,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Eshear'];
units = '$E_{max} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ShearE_max(:,:,:,:));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Evol'];
units = '$E_{vol}=\frac{\delta V}{V} \quad  [\mathrm{mm^3/mm^3}]$';
d=squeeze(Data(series).strain.E_Volumetric(:,:,:,:));
colormap_strain(m,d,suffix,limits,units);

% Strain Lagrange
limits=[-0.7,0.7];
suffix=[foldername,'_Lneg'];
units = '$L_{\lambda_{1}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.L_lambda(:,:,:,:,1));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Lsum'];
units = '$L_{\lambda_{2}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.L_lambda(:,:,:,:,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Lpos'];
units = '$L_{\lambda_{3}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.L_lambda(:,:,:,:,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Lshear'];
units = '$L_{max} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ShearL_max(:,:,:,:));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Lvol'];
units = '$L_{vol}=\frac{\delta V}{V} \quad  [\mathrm{mm^3/mm^3}]$';
d=squeeze(Data(series).strain.L_Volumetric(:,:,:,:));
colormap_strain(m,d,suffix,limits,units);


% Strain Rate
limits=[-1500,1500];
suffix=[foldername,'_SRneg'];
units = '$SR_{\lambda_{1}} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SR_lambda(:,:,:,:,1));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_SRsum'];
units = '$SR_{\lambda_{2}} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SR_lambda(:,:,:,:,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_SRpos'];
units = '$SR_{\lambda_{3}} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SR_lambda(:,:,:,:,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_SRshear'];
units = '$SR_{max} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.ShearSR_max(:,:,:,:));
colormap_strain(m,d,suffix,limits,units);

% Fiber alligned Strain Euler
limits=[-0.4,0.4];
suffix=[foldername,'_Eff'];
units = '$E_{ff} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ED(:,:,:,:,1,1));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Ess'];
units = '$E_{ss} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ED(:,:,:,:,2,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Ett'];
units = '$E_{tt} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ED(:,:,:,:,3,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Efs'];
units = '$E_{fs} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ED(:,:,:,:,1,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Eft'];
units = '$E_{ft} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ED(:,:,:,:,1,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Est'];
units = '$E_{st} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ED(:,:,:,:,2,3));
colormap_strain(m,d,suffix,limits,units);

% Fiber alligned Strain Lagrange
limits=[-0.4,0.4];
suffix=[foldername,'_Lff'];
units = '$L_{ff} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.LD(:,:,:,:,1,1));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Lss'];
units = '$L_{ss} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.LD(:,:,:,:,2,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Ltt'];
units = '$L_{tt} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.LD(:,:,:,:,3,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Lfs'];
units = '$L_{fs} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.LD(:,:,:,:,1,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Lft'];
units = '$L_{ft} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.LD(:,:,:,:,1,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_Lst'];
units = '$L_{st} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.LD(:,:,:,:,2,3));
colormap_strain(m,d,suffix,limits,units);

% Fiber alligned SRex
% Fiber alligned Strain RAte Lagrange
limits=[-750,750];
suffix=[foldername,'_SRff'];
units = '$SR_{ff} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SRD(:,:,:,:,1,1));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_SRss'];
units = '$SR_{ss} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SRD(:,:,:,:,2,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_SRtt'];
units = '$SR_{tt} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SRD(:,:,:,:,3,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_SRfs'];
units = '$SR_{fs} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SRD(:,:,:,:,1,2));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_SRft'];
units = '$SR_{ft} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SRD(:,:,:,:,1,3));
colormap_strain(m,d,suffix,limits,units);

suffix=[foldername,'_SRst'];


units = '$SR_{st} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SRD(:,:,:,:,2,3));
colormap_strain(m,d,suffix,limits,units);

cd ..


end



mkdir('plots')
cd('plots') 

% GM
struct=cat(1,Data(1).strain_GM,Data(2).strain_GM,Data(3).strain_GM);
strain_mvc_plots(struct,'GM')

% SOL
struct=cat(1,Data(1).strain_SOL,Data(2).strain_SOL,Data(3).strain_SOL);
strain_mvc_plots(struct,'SOL')

cd ..

save('Results.mat','Data', '-v7.3')