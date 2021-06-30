function codirected = vectorCOorients(u,v)

% function to check the angle between vetors
% if direction is same (abs(angle) < pi/2) returns true
% if direction is opp. (abs(angle) > pi/2) returns false

CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));

if ThetaInDegrees>90
    codirected=false;
else
    codirected=true;
end

end