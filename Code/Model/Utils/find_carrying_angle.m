function carry = find_carrying_angle
% To find the carrying angle, we place the DSEM model in the initialised
% position, calculate the rotation matrix between humerus and ulna, and
% decompose into Euler angles z - x - y. z is the carrying angle, and at
% the initial position, x and y are both zero.

% Bony landmarks after the model initialisation
IJ        = [     0.00000     0.00000     0.00000]'*0.01; %  (ijplp) 
SC        = [     1.97000     0.93000     1.37000]'*0.01; %  (scplp) 
AC        = [    16.50940     2.66000     7.18000]'*0.01; %  (acplp) 
TS        = [     7.50000    -1.18394    15.59429]'*0.01; %  (tsplp) 
AI        = [    10.16794   -12.63215    15.65449]'*0.01; %  (aiplp) 
GH        = [    17.08274    -1.91903     6.54596]'*0.01; %  (ghjnt) 
E         = [    17.66160   -31.17442     7.03521]'*0.01; %  (hujnt) 
EM        = [    13.68531   -31.06730     6.54596]'*0.01; %  (emplp) 
EL        = [    20.48017   -30.56597     6.54596]'*0.01; %  (elplp) 
WRIST     = [    24.85545   -55.85400     6.55059]'*0.01; %  (rcjnt) 
HAND      = [    28.01019   -65.08909     6.56505]'*0.01; %  (cghan) 
SU        = [    22.17033   -56.40118     6.54595]'*0.01; %  (suplp) 
RD        = [    27.54164   -55.30639     6.54584]'*0.01; %  (srplp) 

Emid = (EM+EL)/2;

% Humerus c.f.
GH_y = (GH-Emid)/norm(GH-Emid);
GH_z1 = cross(EL-EM,GH_y);
GH_z = GH_z1/norm(GH_z1);
GH_x1 = cross(GH_y,GH_z);
GH_x = GH_x1/norm(GH_x1);

%Transformation matrix from g.c.s. to humerus c.s.:
GH_T1 = [GH_x GH_y GH_z GH;0 0 0 1];
GH_T = inv(GH_T1);

%Ulnar c.f.
HU_y = (EM-SU)/norm(EM-SU);
HU_z1 = cross(EL-EM,HU_y);
HU_z = HU_z1/norm(HU_z1);
HU_x1 = cross(HU_y,HU_z);
HU_x = HU_x1/norm(HU_x1);
%Transformation matrix from g.c.s. to ulnar c.s.:
HU_T1 = [HU_x HU_y HU_z E;0 0 0 1];
HU_T = inv(HU_T1);
%Transformation matrix from humerus c.s. to ulnar c.s.:
HU_T_L = GH_T*HU_T1;

[carry,~,~]=rotzxy(HU_T_L(1:3,1:3));