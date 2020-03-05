%% CS_pattern generation
%
%   example call:
%   [K_pattern,CS_pattern,PDF]=CS2GE(106,256,22,2);
%

function [K_pattern,CS_pattern,pdf]=CS2GE(x_res,y_res,nframes,cs_factor)

%% create undersampling pattern
phase_size=x_res;
readout_size=y_res;
frames=nframes;
cs_factor=1/cs_factor;

file_name = strcat('CS_pattern','_x',num2str(1/cs_factor),'.txt');
file_name_mat = strcat('CS','_x',num2str(1/cs_factor),'.mat');

% these are suggested defaults from BART
%cs_factor=.5;
poly_power=20;

[pdf,~]=genPDF([phase_size,1],poly_power,cs_factor,2,0,1);
K_pattern=zeros(frames,phase_size);
N=zeros(frames,1);
b=zeros(nframes,ceil(phase_size*cs_factor));


for i=1:frames

    [mask,~,n] = genSampling(pdf,3,1); 
    N(i)=n;
    a=find(mask)-1;
    a   
        if size(a,2)>ceil(phase_size*cs_factor)
            mask(a(end)) = 0;
            a(end)=[];
        elseif size(a,2)<ceil(phase_size*cs_factor)
            
            for r=0:x_res
                qq=find(a==r);
                if isempty(qq)
                    a=[a, r];
                break    
                end
            end
            
        end
    a=sort(a)    
  
    b(i,:)=a;
    K_pattern(i,:) = mask;

end

b=b';
b=b(:);
b=b';
size(b)
dlmwrite(file_name,b,'-append','delimiter','\n')

CS_pattern=repmat(K_pattern,[1,1,1,readout_size]);
CS_pattern=squeeze(logical(permute(CS_pattern,[4,2,3,1])));
close all

save(file_name_mat)
end
