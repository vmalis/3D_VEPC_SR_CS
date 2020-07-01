%% -----------------------------------------------------------------------   
%                       STRAIN FOR CS VEPC data
%
%  UC San Diego / October 2019 / Vadim Malis
%
%  description: Top level script 4
%
% ------------------------------------------------------------------------
clc
clear all

%% load data
alldata=rdir('**/Results.mat');

%%
load(alldata(1).name);

for i=1:3
Data(i).strain_GM_peak_norm=dev_struct(Data(i).strain_GM_peak,Data(3).strain_GM_peak);
Data(i).strain_SOL_peak_norm=dev_struct(Data(i).strain_SOL_peak,Data(3).strain_SOL_peak);
end

Res=Data;


for i=2:size(alldata,1)
   
    i
    load(alldata(i).name)
    
    for j=1:3
        Data(j).strain_GM_peak_norm=dev_struct(Data(j).strain_GM_peak,Data(3).strain_GM_peak);
        Data(j).strain_SOL_peak_norm=dev_struct(Data(j).strain_SOL_peak,Data(3).strain_SOL_peak);
    end
    
    Res=[Res,Data];
    
end



clearvars -except Res

for i=1:size(Res,2)
    Res(i).MVC_number=str2num(Res(i).MVC{1});
end

%% sort by MVC
Res=nestedSortStruct2(Res,'MVC_number');

%% split by MVC
[~, ia, ~] = unique([Res.MVC_number]);
data30=Res(ia(1):(ia(2)-1));
data40=Res(ia(2):(ia(3)-1));
data60=Res(ia(3):size(Res,2));

%% average of time -----------------------------

%% GM
q=[data30.strain_GM];

fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end

[GM30M,GM30SD] = mean_struct(q);

q=[data40.strain_GM];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[GM40M,GM40SD] = mean_struct(q);

q=[data60.strain_GM];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[GM60M,GM60SD] = mean_struct(q);

GM=[GM30M,GM40M,GM60M];
GM(1).MVC=30;
GM(2).MVC=40;
GM(3).MVC=60;

GM_SD=[GM30SD,GM40SD,GM60SD];
GM_SD(1).MVC=30;
GM_SD(2).MVC=40;
GM_SD(3).MVC=60;



%% SOL
q=[data30.strain_SOL];

fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end

[SOL30M,SOL30SD] = mean_struct(q);

q=[data40.strain_SOL];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[SOL40M,SOL40SD] = mean_struct(q);

q=[data60.strain_SOL];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[SOL60M,SOL60SD] = mean_struct(q);

SOL=[SOL30M,SOL40M,SOL60M];
SOL(1).MVC=30;
SOL(2).MVC=40;
SOL(3).MVC=60;

SOL_SD=[SOL30SD,SOL40SD,SOL60SD];
SOL_SD(1).MVC=30;
SOL_SD(2).MVC=40;
SOL_SD(3).MVC=60;


%% average of peak-----------------------------

%% GM
q=[data30.strain_GM_peak];

fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end

[GM30M_P,GM30SD_P] = mean_struct(q);

q=[data40.strain_GM_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[GM40M_P,GM40SD_P] = mean_struct(q);

q=[data60.strain_GM_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[GM60M_P,GM60SD_P] = mean_struct(q);

GM_P=[GM30M_P,GM40M_P,GM60M_P];
GM_P(1).MVC=30;
GM_P(2).MVC=40;
GM_P(3).MVC=60;

GM_P_SD=[GM30SD_P,GM40SD_P,GM60SD_P];
GM_P_SD(1).MVC=30;
GM_P_SD(2).MVC=40;
GM_P_SD(3).MVC=60;


%% SOL
q=[data30.strain_SOL_peak];

fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end

[SOL30M_P,SOL30SD_P] = mean_struct(q);

q=[data40.strain_SOL_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[SOL40M_P,SOL40SD_P] = mean_struct(q);

q=[data60.strain_SOL_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[SOL60M_P,SOL60SD_P] = mean_struct(q);

SOL_P=[SOL30M_P,SOL40M_P,SOL60M_P];
SOL_P(1).MVC=30;
SOL_P(2).MVC=40;
SOL_P(3).MVC=60;

SOL_P_SD=[SOL30SD_P,SOL40SD_P,SOL60SD_P];
SOL_P_SD(1).MVC=30;
SOL_P_SD(2).MVC=40;
SOL_P_SD(3).MVC=60;





%% average of normalized peak-----------------------------

%% GM
q=[data30.strain_GM_peak_norm];

fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end

[GM30M_PN,GM30SD_PN] = mean_struct(q);

q=[data40.strain_GM_peak_norm];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[GM40M_PN,GM40SD_PN] = mean_struct(q);

q=[data60.strain_GM_peak_norm];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[GM60M_PN,GM60SD_PN] = mean_struct(q);

GM_PN=[GM30M_PN,GM40M_PN,GM60M_PN];
GM_PN(1).MVC=30;
GM_PN(2).MVC=40;
GM_PN(3).MVC=60;

GM_PN_SD=[GM30SD_PN,GM40SD_PN,GM60SD_PN];
GM_PN_SD(1).MVC=30;
GM_PN_SD(2).MVC=40;
GM_PN_SD(3).MVC=60;


%% SOL
q=[data30.strain_SOL_peak_norm];

fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end

[SOL30M_PN,SOL30SD_PN] = mean_struct(q);

q=[data40.strain_SOL_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[SOL40M_PN,SOL40SD_PN] = mean_struct(q);

q=[data60.strain_SOL_peak_norm];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
[SOL60M_PN,SOL60SD_PN] = mean_struct(q);

SOL_PN=[SOL30M_PN,SOL40M_PN,SOL60M_PN];
SOL_PN(1).MVC=30;
SOL_PN(2).MVC=40;
SOL_PN(3).MVC=60;

SOL_PN_SD=[SOL30SD_PN,SOL40SD_PN,SOL60SD_PN];
SOL_PN_SD(1).MVC=30;
SOL_PN_SD(2).MVC=40;
SOL_PN_SD(3).MVC=60;


clearvars -except SOL_P_SD SOL_PN_SD SOL_P SOL_PN SOL SOL_SD GM_P_SD GM_P GM GM_SD GM_PN_SD GM_PN GM GM_SD data30 data40 data60

%%make struct for plots
gm  = meanANDsd(GM,GM_SD);
sol = meanANDsd(SOL,SOL_SD);


%% plot time average
strain_mvc_plots(nestedSortStruct(gm,'MVC',-1),'GM')
strain_mvc_plots(nestedSortStruct(sol,'MVC',-1),'SOL')
clearvars('SOL', 'SOL_SD', 'GM', 'GM_SD','gm','sol');









%% do stats 

%%% prepare for anova


%------------absolute SOL----------------
q=[data30.strain_SOL_peak];
fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
SOL30P=q;
for i=1:size(SOL30P,2)
SOL30P(i).MVC=30;
end

q=[data40.strain_SOL_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
SOL40P=q;
for i=1:size(SOL40P,2)
SOL40P(i).MVC=40;
end

q=[data60.strain_SOL_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
SOL60P=q;
for i=1:size(SOL60P,2)
SOL60P(i).MVC=60;
end


%------------normalized SOL----------------
q=[data30.strain_SOL_peak_norm];
fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
SOL30PN=q;
for i=1:size(SOL30PN,2)
SOL30PN(i).MVC=30;
end

q=[data40.strain_SOL_peak_norm];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
SOL40PN=q;
for i=1:size(SOL40PN,2)
SOL40PN(i).MVC=40;
end

q=[data60.strain_SOL_peak_norm];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
SOL60PN=q;
for i=1:size(SOL60PN,2)
SOL60PN(i).MVC=60;
end




%------------absolute GM----------------
q=[data30.strain_GM_peak];
fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
GM30P=q;
for i=1:size(GM30P,2)
GM30P(i).MVC=30;
end

q=[data40.strain_GM_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
GM40P=q;
for i=1:size(GM40P,2)
GM40P(i).MVC=40;
end

q=[data60.strain_GM_peak];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
GM60P=q;
for i=1:size(GM60P,2)
GM60P(i).MVC=60;
end

%------------normalized GM----------------

q=[data30.strain_GM_peak_norm];
fnames=fieldnames(q);
sd_fields = find(contains(fnames,'sd'));

for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
GM30PN=q;
for i=1:size(GM30PN,2)
GM30PN(i).MVC=30;
end

q=[data40.strain_GM_peak_norm];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
GM40PN=q;
for i=1:size(GM40PN,2)
GM40PN(i).MVC=40;
end

q=[data60.strain_GM_peak_norm];
for i=1:size(sd_fields,1)
    q=rmfield(q,fnames{sd_fields(i)});
end
GM60PN=q;
for i=1:size(GM60PN,2)
GM60PN(i).MVC=60;
end




GMP=[GM60P,GM40P,GM30P];
GMPN=[GM60PN,GM40PN,GM30PN];
SOLP=[SOL60P,SOL40P,SOL30P];
SOLPN=[SOL60PN,SOL40PN,SOL30PN];


clearvars("SOL60P","SOL40P","SOL30P","GM60P","GM40P","GM30P","SOL60PN","SOL40PN","SOL30PN","GM60PN","GM40PN","GM30PN","fnames","q","i","data30","data40","data60","sd_fields")

% last step before stats
% flattenen the structure so all fields are size 1, if necessary new fields
% are created


GM_peak  = flattenSRstruct(GMP);                %ANOVA READY Peak data
SOL_peak = flattenSRstruct(SOLP);               %ANOVA READY Peak data
GM_peak_norm  = flattenSRstruct(GMPN);          %ANOVA READY Peak data
SOL_peak_norm = flattenSRstruct(SOLPN);         %ANOVA READY Peak data


%%
GM       = flattenSRstruct(GM_P);               %MEAN TABLE
GM_SD    = flattenSRstruct(GM_P_SD);            %SD TABLE
GM_N       = flattenSRstruct(GM_PN);            %MEAN TABLE
GM_N_SD    = flattenSRstruct(GM_PN_SD);         %SD TABLE


SOL       = flattenSRstruct(SOL_P);             %MEAN TABLE
SOL_SD    = flattenSRstruct(SOL_P_SD);          %SD TABLE
SOL_N       = flattenSRstruct(SOL_PN);          %MEAN TABLE
SOL_N_SD    = flattenSRstruct(SOL_PN_SD);       %SD TABLE

%get rid of evrythin except these structures
clearvars -except SOL_peak SOL_peak_norm GM_peak GM_peak_norm GM GM_SD SOL SOL_SD GM_N GM_N_SD SOL_N SOL_N_SD
 

%% ANOVA
% in here just grab a structure/structure and do anova
%
% % the results are represented as following

GM_ANOVA  = anovaCSStrain(GM_peak);
SOL_ANOVA = anovaCSStrain(SOL_peak);







