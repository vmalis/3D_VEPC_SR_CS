%_____________________________________________________
% Subroutine to perform SR to DTI rotation
%_____________________________________________________
% written by Vadim Malis
% 10/17 at UCSD RIL

function [SR_D]=dti_rotation(rotation_matrix,tensor_to_rotate)

SR_D=zeros(size(tensor_to_rotate));


ii=size(tensor_to_rotate,1);
jj=size(tensor_to_rotate,2);
kk=size(tensor_to_rotate,3);
tt=size(tensor_to_rotate,4);

    for i=1:ii
        for j=1:jj
            for k=1:kk
                for t=1:tt
                                        
                    R = squeeze(rotation_matrix(i,j,k,:,:));
                    Q = flip(squeeze(tensor_to_rotate(i,j,k,t,:,:)),2);
                    SR_D(i,j,k,t,:,:)=R*Q*R';

                     end
                    
                end
            end
        end
    end
    
