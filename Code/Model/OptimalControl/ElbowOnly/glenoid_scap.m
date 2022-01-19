function Rginv = glenoid_scap
% finds rotation matrix between scapular coordinate frame and glenoid
% orientation based on the DSEM cadaver file
% all positions are in the global coordinate frame, in cm
%
% Dimitra Blana, March 2012

% Get positions from dsp file:
% position glenohumeral joint
GH_centre = [17.08   -1.92    6.55];
% mid-point articular surface glenoid
glenoid_centre=[15.14   -1.98    7.81];
% In vivo palpated bony landmark at the scapula (AA):
AA = [18.26    0.75   10.56];
% In vivo palpated bony landmark at the scapula (TS):
TS = [7.50   -1.17   15.60];
% In vivo palpated bony landmark at the scapula (AI):
AI=[10.16  -12.62   15.67];

% Find scapular coordinate frame:
% local x-axis : AA to TS               										
% local z-axis : perpendicular to x and the line connecting AA and AI     	
% local y-axis : perpendicular to z and x				                       	
Xs = (AA-TS) / norm(AA-TS);
Zs = cross(Xs,(AA-AI)); 
Zs = Zs/norm(Zs);
Ys = cross(Zs,Xs);
S = [Xs;Ys;Zs];

%% Find vector from glenoid to GH centre in the global frame:
glen_scap_v = glenoid_centre - GH_centre;
% in scapular frame:
glen_scap = S*glen_scap_v';

% find polar angles:
thetar = asin(-glen_scap(2));
if ~(sqrt(glen_scap(1)^2+glen_scap(3)^2))
    phir = 0.0;
else
    phir = asin(glen_scap(3)/(sqrt(glen_scap(1)^2+glen_scap(3)^2)));
end

% calculate orientation matrix of glenoid, Rg, and inverse
Rg=roty(phir)*rotz(thetar);
Rginv = Rg';
    

function [Ry]=roty(th)
% creates rotation matrix
% for rotations or th radians around the y axis
Ry(1,1)=cos(th);
Ry(1,3)=sin(th);
Ry(2,2)=1;
Ry(3,1)=-sin(th);
Ry(3,3)= cos(th);

function [Rz]=rotz(th)
% calculates the rotation matrix
% for rotations of th radians around the z axis

Rz(1,1)=cos(th);
Rz(1,2)=-sin(th);
Rz(2,1)= sin(th);
Rz(2,2)= cos(th);
Rz(3,3)=1;
