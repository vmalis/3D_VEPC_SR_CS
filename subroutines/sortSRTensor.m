function B = sortSRTensor(A)

% function to sort secondary a tertiary elements of tensor

% input is TENSOR
%
%               | 1,1   2,1    3,1 |
%     TENSOR =  | 1,2   2,2    3,2 |
%               | 1,3   2,3    3,3 |
%
%
%     1,1  is always ff
%               
%     if abs(2,2) > abs (3,3) we have to switch indecies 2 and 3
%

B=A;

    if abs(A(2,2))>abs(A(3,3))

        B(2,2) = A(3,3);
        B(3,3) = A(2,2);

        B(2,1) = A(3,1);
        B(1,2) = B(2,1);

        B(3,1) = A(2,1);
        B(1,3) = B(3,1);

        % elemnts 3,2 and 2,3 are same so no need to swap

    end


end