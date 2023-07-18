function [Tthorax,Tclavicle,Tscapula,Thumerus,Tulna,Tradius,Thand,humcorr] = cfb(option)
% definition of coordinate frames
% if option = 1, the output is the transformation matrices from the global
% to the local coordinate frames
% if option = 2, the output is the transformation matrices from the child
% segment to the parent segment coordinate frame
% if option = 3, the output is the transformation matrices from the local
% to the global coordinate frames
% if option = 4, the output is the transformation matrices from the
% parent to the child coordinate frames

% bony landmarks
IJ = [0;0;0];
SC = [0.0197;0.0093;0.0137];
AC = [0.1651;0.0266;0.0718];
AA = [0.1826;0.0075;0.1056];
TS = [0.0750;-0.0117;0.1560];
AI = [0.1016;-0.1262;0.1567];
HH = [0.1686;-0.0175;0.0651]; %center of humeral head
GH = [0.1708;-0.0192;0.0655]; %glenohumeral joint


% the humeral head is dislocated: this is the translation all points on
% humerus and beyond have to be corrected by
humcorr = HH-GH;

EM = [0.1566;-0.3079;0.1056]-humcorr;       %epic. mid
EL = [0.2153;-0.3021;0.0715]-humcorr;       %epic. lat.
E = (EM+EL)/2;                              %midpoint
OL = [0.1897;-0.3023;0.1085]-humcorr;       %olecranon
SU = [0.2154;-0.5571;0.0339]-humcorr;       %styloid process at the ulna
SR = [0.1729;-0.5501;0.0001]-humcorr;       %styloid process at the radius
WR = (SU+SR)/2;
CR = [0.1822; -0.6473; -0.0075]-humcorr;    %position of force application on the hand

% joint positions
SC_origin = [0.0014; -0.0152; 0.0028];
AC_origin = AA;
GH_origin = [0.1708; -0.0192; 0.0655];
HU_origin = [0.1936; -0.3080; 0.0902]-humcorr;
UR_origin = [0.2011; -0.3227; 0.0687]-humcorr;
RC_origin = [0.1942; -0.5536; 0.0170]-humcorr;

%Thorax c.f.: alligned with the global c.f.
TH_x = [1;0;0];
TH_y = [0;1;0];
TH_z = [0;0;1];
TH_origin = IJ; 
TH_T = [TH_x TH_y TH_z TH_origin;0 0 0 1];

%Clavicular c.f.
SC_x = (AC-SC)/norm(AC-SC);
SC_z1 = cross(SC_x,TH_y);
SC_z = SC_z1/norm(SC_z1);
SC_y1 = cross(SC_z,SC_x);
SC_y = SC_y1/norm(SC_y1);

%Transformation matrix from g.c.s. to clavicular c.s.:
SC_T1 = [SC_x SC_y SC_z SC_origin;0 0 0 1];
SC_T = inv(SC_T1);
%(P at c.c.s. = SC_T * P at g.c.s.)
% Transformation matrix from clavicular c.s. to thorax c.s.:
SC_T_L = TH_T*SC_T1;

%Scapular c.f.
AC_x = (AA-TS)/norm(AA-TS);
AC_z1 = cross(AI-AA,AC_x);
AC_z = AC_z1/norm(AC_z1);
AC_y1 = cross(AC_z,AC_x);
AC_y = AC_y1/norm(AC_y1);

%Transformation matrix from g.c.s. to scapular c.s.:
AC_T1 = [AC_x AC_y AC_z AC_origin;0 0 0 1];
AC_T = inv(AC_T1);
%(P at s.c.s. = AC_T * P at g.c.s.)
%Transformation matrix from clavicular c.s. to scapular c.s.:
AC_T_L = SC_T*AC_T1;

% Humerus c.f.
GH_y = (GH-E)/norm(GH-E);
GH_z1 = cross(EL-EM,GH_y);
GH_z = GH_z1/norm(GH_z1);
GH_x1 = cross(GH_y,GH_z);
GH_x = GH_x1/norm(GH_x1);

%Transformation matrix from g.c.s. to humerus c.s.:
GH_T1 = [GH_x GH_y GH_z GH_origin;0 0 0 1];
GH_T = inv(GH_T1);
%Transformation matrix from scapular c.s. to humerus c.s.:
GH_T_L = AC_T*GH_T1;

%Ulnar c.f.
HU_y = (EM-SU)/norm(EM-SU);
HU_z1 = cross(EL-EM,HU_y);
HU_z = HU_z1/norm(HU_z1);
HU_x1 = cross(HU_y,HU_z);
HU_x = HU_x1/norm(HU_x1);
%Transformation matrix from g.c.s. to ulnar c.s.:
HU_T1 = [HU_x HU_y HU_z HU_origin;0 0 0 1];
HU_T = inv(HU_T1);
%Transformation matrix from humerus c.s. to ulnar c.s.:
HU_T_L = GH_T*HU_T1;

% Radius c.f.
UR_y = (EL-SR)/norm(EL-SR);
UR_z1 = cross(SR-SU,UR_y);
UR_z = UR_z1/norm(UR_z1);
UR_x1 = cross(UR_y,UR_z);
UR_x = UR_x1/norm(UR_x1);
%Transformation matrix from g.c.s. to radius c.s.:
UR_T1 = [UR_x UR_y UR_z UR_origin;0 0 0 1];
%UR_T1 = [UR_x UR_y UR_z SR;0 0 0 1];
UR_T = inv(UR_T1);
%Transformation matrix from ulnar c.s. to radius c.s.:
UR_T_L = HU_T*UR_T1;

% Hand c.f.: alligned with the radius c.f.
%Transformation matrix from g.c.s. to hand c.s.:
RC_T1 = [UR_x UR_y UR_z RC_origin;0 0 0 1];
RC_T = inv(RC_T1);
%Transformation matrix from radius c.s. to hand c.s.:
RC_T_L = UR_T*RC_T1;


switch option
    case 1
        Tthorax = TH_T;
        Tclavicle = SC_T;
        Tscapula = AC_T;
        Thumerus = GH_T;
        Tulna = HU_T;
        Tradius = UR_T;
        Thand = RC_T;
    case 2
        Tthorax = TH_T;
        Tclavicle = SC_T_L;
        Tscapula = AC_T_L;
        Thumerus = GH_T_L;
        Tulna = HU_T_L;
        Tradius = UR_T_L;
        Thand = RC_T_L;
    case 3
        Tthorax = TH_T;
        Tclavicle = SC_T1;
        Tscapula = AC_T1;
        Thumerus = GH_T1;
        Tulna = HU_T1;
        Tradius = UR_T1;
        Thand = RC_T1;
    case 4
        Tthorax = TH_T;
        Tclavicle = inv(SC_T_L);
        Tscapula = inv(AC_T_L);
        Thumerus = inv(GH_T_L);
        Tulna = inv(HU_T_L);
        Tradius = inv(UR_T_L);
        Thand = inv(RC_T_L);
end

