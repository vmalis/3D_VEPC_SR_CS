%% Univariate analysis

%% read data
data=csv2struct("CS_MVC.xlsx");
fnames = {'MVCpt', 'trueMVCpt', 'Force'};
fnamesRx=fieldnames(data);
fnamesRx(1:7)=[];

%create combinations

%different combinations
DATA(1).name="All";
DATA(1).input=data;

%young GM
idx = (ismember({data.Age}, {'O'})+ ismember({data.Muscle}, {'SOL'}));
idx(idx>1)=1;
DATA(2).name="Young GM";
DATA(2).input= data(~idx);

%senior GM
idx = (ismember({data.Age}, {'Y'})+ ismember({data.Muscle}, {'SOL'}));
idx(idx>1)=1;
DATA(3).name="Senior GM";
DATA(3).input= data(~idx);

%young SOL
idx = (ismember({data.Age}, {'O'})+ ismember({data.Muscle}, {'GM'}));
idx(idx>1)=1;
DATA(4).name="Young SOL";
DATA(4).input= data(~idx);

%senior SOL
idx = (ismember({data.Age}, {'Y'})+ ismember({data.Muscle}, {'GM'}));
idx(idx>1)=1;
DATA(5).name="Senior SOL";
DATA(5).input = data(~idx);

%young
idx = ismember({data.Age}, {'O'});
DATA(6).name="Young";
DATA(6).input = data(~idx);

%senior
idx = ismember({data.Age}, {'Y'});
DATA(7).name="Senior";
DATA(7).input = data(~idx);


for i = 1:3 %for dependet
    for j = 1:7 %for different combinations
        DATA(j).statsMVCpt = CSMVC_univariate([DATA(j).input.(fnames{1})],DATA(j).input,fnamesRx);
        DATA(j).statstrueMVCpt = CSMVC_univariate([DATA(j).input.(fnames{2})],DATA(j).input,fnamesRx);
        DATA(j).statsForce = CSMVC_univariate([DATA(j).input.(fnames{3})],DATA(j).input,fnamesRx);
    end
end



