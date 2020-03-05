%% Diffusion preprocessing for 3D Strain:
%
% all in one script for Diffusion data processing
%
% _____________________________________________________
% written by Vadim Malis
% 07/19 at UCSD RIL



%% run eddy current and field map

[DWI,Baseline,Grad,Location,info]=fsl_vepc;

DWI=double(DWI);
Baseline=double(Baseline);
% geometry amd b value information for merged volume
info.z=[Location(end) Location(1)];
info.nx=size(DWI,2);
info.ny=size(DWI,1);
info.nz=size(DWI,3);

clearvars -except Baseline DWI Grad Location info 

%% perform LMMSE filtering
%get sigma for noise estimation in the baseline
sigma=RicianSTD(squeeze(DWI(:,:,:,1)));
%filter
DWI=jaLMMSEDWI(DWI,Grad,sigma*5);


%% estimate tensor
[DTI,Lambda,FA,ADC,Vector,CP]=diffusion_tensor(DWI,Baseline,Grad,info.b);


%% saving eigenvector colormaps as a png for preview
EV1=abs(permute(squeeze(Vector(:,:,:,:,3)),[1,2,4,3]));
EV1=EV1(:,:,[1,3,2],:);
EV2=abs(permute(squeeze(Vector(:,:,:,:,2)),[1,2,4,3]));
EV2=EV2(:,:,[1,3,2],:);
EV3=abs(permute(squeeze(Vector(:,:,:,:,1)),[1,2,4,3]));
EV3=EV3(:,:,[1,3,2],:);

montage(EV1)
export_fig('ev1.png', '-png','-m16')
close

montage(EV2)
export_fig('ev2.png', '-png','-m16')
close

montage(EV3)
export_fig('ev3.png', '-png','-m16')
close


% colormaps for ADC, FA, CP

cmap = jet(201);

dti_montage('ADC',ADC,cmap,1,[0 0.003]);
dti_montage('FA',FA,cmap,1,[0 .5]);
dti_montage('CP',CP,cmap,1,[0 .5]);


DTI_data.ID=info.ID;
DTI_data.header=info;
DTI_data.location=Location;
DTI_data.DTI_eigenvector=Vector;
DTI_data.DTI_eigenvalue=Lambda;
DTI_data.DTI=DTI;
DTI_data.FA=FA;
DTI_data.ADC=ADC;
DTI_data.CP=CP;

save('DTI.mat', 'DTI_data');



close all
clear all
clc
cd ..
cd ..



