% Script to convert DSEM data to Opensim. This focuses on the humerus and
% forearm, as the shoulder data has already been checked.
%
% The script calls cfb.m, which calculates transformation matrices between
% coordinates frames based on bony landmarks, and find_carrying_angle.m
% that calculates the forearm carrying angle. 
%
% It loads DSEM_muscles.mat, which contains the muscle attachment points
% from the DSEM model file l1091.dsp. (To do: find the code that generated
% this .mat file, as coracobrachialis appears to be missing.)
% 
% The script requires rotxyz.m and rotzxy.m (for the carrying angle).

%% Set up output file

outfile = fopen('DSEM_to_OpenSim.txt', 'w');
fprintf(outfile,'This file contains data from the DSEM model\n');
fprintf(outfile,'converted to the format required for Opensim\n\n');

%% Get carrying angle
carry = find_carrying_angle;
disp_str = ['The carrying angle is ' num2str(carry) ' radians (' num2str(carry*180/pi) ' degrees)'];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);

%% Tranformation matrix from the child segment to the parent segment coordinate frame
[Tthorax,Tclavicle,Tscapula,Thumerus,Tulna,Tradius,Thand,~] = cfb(2);

% Get location of all joint centres in parent frames
SCjoint = Tclavicle(1:3,4)';
ACjoint = Tscapula(1:3,4)';
GHjoint = Thumerus(1:3,4)';
HUjoint = Tulna(1:3,4)';
URjoing = Tradius(1:3,4)';
RCjoing = Thand(1:3,4)';

%% Get locations of bony landmarks in local frames
% Tranformation matrix from the global to the local coordinate frames
[Tthorax,Tclavicle,Tscapula,Thumerus,Tulna,Tradius,Thand,humcorr] = cfb(1);

SC = [1.97    0.93    1.37]'/100;
AA = [18.26    0.75   10.56]'/100;
AI = [10.16  -12.62   15.67]'/100;
TS = [7.50   -1.17   15.60]'/100;
AC = [16.51    2.66    7.18]'/100;
GH = [16.86   -1.75    6.51]'/100 - humcorr;
EM = [15.66  -30.79   10.56]'/100 - humcorr;
EL = [21.53  -30.21    7.15]'/100 - humcorr;
US = [21.54  -55.71    3.39]'/100 - humcorr;
RS = [17.29  -55.01    0.01]'/100 - humcorr;

SC_local = Tclavicle*[SC;1];
AA_local = Tscapula*[AA;1];
AI_local = Tscapula*[AI;1];
TS_local = Tscapula*[TS;1];
AC_local = Tscapula*[AC;1];
GH_local = Thumerus*[GH;1];
EM_local = Thumerus*[EM;1];
EL_local = Thumerus*[EL;1];
US_local = Tulna*[US;1];
RS_local = Tradius*[RS;1];

%% Elbow flexion-extension and pronation-supination axes

% Tranformation matrix from the parent segment to the child segment coordinate frame
% [~,~,~,~,Tulna1,Tradius1,~,~] = cfb(4);
% 
% flex_ext_axis = [0.8531   0.0183  -0.5108]';
% flexext = Tulna1*[flex_ext_axis;1];
% flex_ext_axis_local = flexext(1:3)/norm(flexext(1:3));
% disp(['Elbow flexion-extension axis: ' num2str(flex_ext_axis_local')]);

% Set the flexion-extension axis as the line from EL to EM
EL_in_ulna = [0.0280109 -0.00501545 -0.0049];
EM_in_ulna = [-0.0379636 0.0118005 -0.0049];
flex_ext_axis_local = (EL_in_ulna - EM_in_ulna)/norm(EL_in_ulna-EM_in_ulna);
disp_str = ['Elbow flexion-extension axis: ' num2str(flex_ext_axis_local)];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

% pro_sup_axis = [-0.0604   0.9874   0.1465]';
% prosup = Tradius1*[pro_sup_axis;1];
% pro_sup_axis_local = prosup(1:3)/norm(prosup(1:3));
% disp(['Pronation-supination axis: ' num2str(pro_sup_axis_local')]);

% Set the pronation-supination axis as the line from ulnar styloid to
% capitulum humeri
CH_in_radius = [-0.00287672 0.0155471 -0.00434069]; % New point Capitulum Humeri estimated in radius frame
US_in_radius = [-0.0488 -0.2323 0.007]; 
EL_in_radius = [0.0169109 0.0161846 0.0078];
pro_sup_axis_local = (CH_in_radius - US_in_radius)/norm(CH_in_radius - US_in_radius);
disp_str = ['Pronation-supination axis: ' num2str(pro_sup_axis_local)];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);


%% all muscle points 

mus = load('DSEM_muscles');

for imus=1:length(mus.muscles)
    for ielem=1:length(mus.muscles{imus}.elem)
        for ipt=1:size(mus.muscles{imus}.elem{ielem}.points,1)
            coords = mus.muscles{imus}.elem{ielem}.points(ipt,2:4);     
            if (norm(coords - [0.0999 0.0999 0.0999])<0.01 || norm(coords - [0.1021 0.0982 0.1003])<0.01), continue; end
            switch mus.muscles{imus}.elem{ielem}.points(ipt,1)
                case 24 % scapula
                    coords_loc = (Tscapula*([coords,1]'))';
                    disp_str = [mus.muscles{imus}.name num2str(ielem) ' point ' num2str(ipt) ' on scapula:' num2str(coords_loc(1:3))];
                    disp(disp_str);
                    fprintf(outfile,'%s\n',disp_str);
                case 25 % humerus
                    coords_loc = (Thumerus*([coords,1]'))';
                    disp_str = [mus.muscles{imus}.name num2str(ielem) ' point ' num2str(ipt) ' on humerus:' num2str(coords_loc(1:3))];
                    disp(disp_str);
                    fprintf(outfile,'%s\n',disp_str);
                case 26 % ulna
                    coords_loc = (Tulna*([coords,1]'))';
                    disp_str = [mus.muscles{imus}.name num2str(ielem) ' point ' num2str(ipt) ' on ulna:' num2str(coords_loc(1:3))];
                    disp(disp_str);
                    fprintf(outfile,'%s\n',disp_str);
                case 27 % radius
                    coords_loc = (Tradius*([coords,1]'))';
                    disp_str = [mus.muscles{imus}.name num2str(ielem) ' point ' num2str(ipt) ' on radius:' num2str(coords_loc(1:3))];
                    disp(disp_str);
                    fprintf(outfile,'%s\n',disp_str);
            end
        end
    end
end

% Coracobrachialis is missing...? Just humerus attachment of second element for now:
coords = ([17.20  -18.35    8.10]'/100 - humcorr)';
coords_loc = (Thumerus*([coords,1]'))';
disp_str = ['Coracobrachialis_2 point 2 on humerus:' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

%% wrapping objects

fprintf(outfile,'\nWrapping Objects:\n');

% REM ELLIPSO surfacenr mx my mz ax ay az pp po
% REM mx, my, mz: coordinates of centre of ellipsoid
% REM ax, ay, az: lengths of the axes of the ellipsoid
% REM pp, po: position and orientation nodes
% 
% ELLIPSO  1   0.00 -14.86   5.91  14.70  20.79   9.44 23  4
% 
% REM BALL surfacenr mx my mz r pp po
% REM mx, my, mz: coordinates of centre of ball
% REM r: radius of the ball
% REM pp, po: position and orientation nodes
% 
% BALL  2  17.19  -2.03   6.29   2.72 26 13
% BALL  3  16.86  -1.75   6.51   2.31 26 13

coords = [17.19  -2.03   6.29]'/100 - humcorr;
coords_loc = (Thumerus*([coords;1]))';
disp_str = [' Centre of ball "deltoid" on humerus: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

coords = [16.86  -1.75   6.51]'/100 - humcorr;
coords_loc = (Thumerus*([coords;1]))';
disp_str = [' Centre of ball "inf_sup_sub" on humerus: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);


% REM CYLINDER surfacenr dx dy dz sx sy sz r pp po
% REM dx, dy, dz: unit vectors in direction of axes
% REM sx, sy, sz: coordinates of a point on the long axis
% REM r: radius of the cylinder
% REM pp, po: position and orientation nodes
% 
% CYLINDER  4  0.0602 -0.9949  0.0810  17.75  -9.22   6.49  1.00 26 13

vec = [0.0602; -0.9949;  0.0810];
v_local = Thumerus*([vec;1]);
z = v_local(1:3)/norm(v_local(1:3));
x = cross(z,[1;0;0]);
y = cross(z,x);
x = cross(y,z);
R = [x,y,z];
[angles(1), angles(2), angles(3)] = rotxyz(R);
disp_str = [' Rotation of cylinder "infra_cyl" on humerus: ' num2str(angles)];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

coords = [17.75  -9.22   6.49]'/100 - humcorr;
coords_loc = (Thumerus*([coords;1]))';
disp_str = [' Translation of cylinder "infra_cyl" on humerus: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);

vec = [0.0993;  0.9003;  0.4239];
v_local = Thumerus*([vec;1]);
z = v_local(1:3)/norm(v_local(1:3));
x = cross(z,[1;0;0]);
y = cross(z,x);
x = cross(y,z);
R = [x,y,z];
[angles(1), angles(2), angles(3)] = rotxyz(R);
disp_str = [' Rotation of cylinder "ter_lat_pec" on humerus: ' num2str(angles)];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

coords = [19.92 -35.32   5.79]'/100 - humcorr;
coords_loc = (Thumerus*([coords;1]))';
disp_str = [' Translation of cylinder "ter_lat_pec" on humerus: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);


% ulna

% CYLINDER  7 -0.8531 -0.0183  0.5108  19.36 -30.80   9.02  1.50 27 16

vec = [-0.8531; -0.0183; 0.5108];
v_local = Tulna*([vec;1]);
z = v_local(1:3)/norm(v_local(1:3));
x = cross(z,[1;0;0]);
y = cross(z,x);
x = cross(y,z);
R = [x,y,z];
[angles(1), angles(2), angles(3)] = rotxyz(R);
disp_str = [' Rotation of cylinder "elbow" on ulna: ' num2str(angles)];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

coords = [19.36 -30.80   9.02]'/100 - humcorr;
coords_loc = (Tulna*([coords;1]))';
disp_str = [' Translation of cylinder "elbow" on ulna: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);

% CYLINDER  6  0.8531  0.0183 -0.5108  19.36 -30.80   9.02  1.90 27 16

vec = [0.8531;  0.0183; -0.5108];
v_local = Tulna*([vec;1]);
z = v_local(1:3)/norm(v_local(1:3));
x = cross(z,[1;0;0]);
y = cross(z,x);
x = cross(y,z);
R = [x,y,z];
[angles(1), angles(2), angles(3)] = rotxyz(R);
disp_str = [' Rotation of cylinder "elbow_circ" on ulna: ' num2str(angles)];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

coords = [19.36 -30.80   9.02]'/100 - humcorr;
coords_loc = (Tulna*([coords;1]))';
disp_str = [' Translation of cylinder "elbow_circ" on ulna: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);

z = flex_ext_axis_local';
x = cross(z,[1;0;0]);
y = cross(z,x);
x = cross(y,z);
R = [x,y,z];
[angles(1), angles(2), angles(3)] = rotxyz(R);
disp_str = [' Rotation of cylinder "elbow_circ" on ulna (aligned with flex/ext axis: ' num2str(angles)];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);

% CYLINDER  9 -0.1157  0.9547  0.2740  20.77 -51.88   3.88  0.70 27 16

vec = [-0.1157; 0.9547; 0.2740];
v_local = Tulna*([vec;1]);
z = v_local(1:3)/norm(v_local(1:3));
x = cross(z,[1;0;0]);
y = cross(z,x);
x = cross(y,z);
R = [x,y,z];
[angles(1), angles(2), angles(3)] = rotxyz(R);
disp_str = [' Rotation of cylinder "pron_quad" on ulna: ' num2str(angles)];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

coords = [20.77 -51.88   3.88]'/100 - humcorr;
coords_loc = (Tulna*([coords;1]))';
disp_str = [' Translation of cylinder "pron_quad" on ulna: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);


% radius

% CYLINDER 10 -0.1481 -0.8839 -0.4436  19.28 -40.60   3.49  0.71 28 19

% vec = [-0.1481; -0.8839; -0.4436];
% v_local = Tradius*([vec;1]);
% z = v_local(1:3)/norm(v_local(1:3));
% x = cross(z,[1;0;0]);
% y = cross(z,x);
% x = cross(y,z);
% R = [x,y,z];
% [angles(1), angles(2), angles(3)] = rotxyz(R);
% disp_str = [' Rotation of cylinder "ulna_radius" on radius: ' num2str(angles)];
% disp(disp_str);
% fprintf(outfile,'%s\n',disp_str);
% 
% coords = [19.28 -40.60   3.49]'/100 - humcorr;
% coords_loc = (Tradius*([coords;1]))';
% disp_str = [' Translation of cylinder "ulna_radius" on radius: ' num2str(coords_loc(1:3))];
% disp(disp_str);
% fprintf(outfile,'%s\n\n',disp_str);

% CYLINDER  5  0.0993  0.9003  0.4239  19.92 -35.32   5.79  0.92 28 19

vec = [0.0993;  0.9003; 0.4239];
v_local = Tradius*([vec;1]);
z = v_local(1:3)/norm(v_local(1:3));

z = pro_sup_axis_local';
x = cross(z,[1;0;0]);
y = cross(z,x);
x = cross(y,z);
R = [x,y,z];
[angles(1), angles(2), angles(3)] = rotxyz(R);
disp_str = [' Rotation of cylinder "bic_radius" on radius (aligned with pro/sup axis): ' num2str(angles)];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);


coords = [19.92 -35.32   5.79]'/100 - humcorr;
coords_loc = (Tradius*([coords;1]))';
disp_str = [' Translation of cylinder "bic_radius" on radius: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);

% Moved out slightly:
coords_loc = [-0.004 -0.061 0.012];
disp_str = [' Translation of cylinder "bic_radius" on radius (moved slightly): ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n\n',disp_str);

% CYLINDER  8  0.0186  0.9767  0.2136  19.79 -43.87   3.87  0.90 28 19

% vec = [0.0186;  0.9767;  0.2136];
% v_local = Tradius*([vec;1]);
% z = v_local(1:3)/norm(v_local(1:3));
% x = cross(z,[1;0;0]);
% y = cross(z,x);
% x = cross(y,z);
% R = [x,y,z];
% [angles(1), angles(2), angles(3)] = rotxyz(R);
% disp_str = [' Rotation of cylinder "supinator" on radius: ' num2str(angles)];
% disp(disp_str);
% fprintf(outfile,'%s\n',disp_str);

supinator_axis = (EL_in_radius' - RS_local(1:3))/norm(EL_in_radius' - RS_local(1:3));
v_local = supinator_axis;
z = v_local(1:3)/norm(v_local(1:3));
x = cross(z,[1;0;0]);
y = cross(z,x);
x = cross(y,z);
R = [x,y,z];
[angles(1), angles(2), angles(3)] = rotxyz(R);
disp_str = [' Rotation of cylinder "supinator" on radius: ' num2str(angles)];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

coords = [19.79 -43.87   3.87]'/100 - humcorr;
coords_loc = (Tradius*([coords;1]))';
disp_str = [' Translation of cylinder "supinator" on radius: ' num2str(coords_loc(1:3))];
disp(disp_str);
fprintf(outfile,'%s\n',disp_str);

fclose(outfile);

% Pronator quadratus is too long...? (has moment arm that is too large)
att1 = [-0.048 -0.213 -0.005];

old_att2 = [-0.015615 -0.21256 0.0087301];
new_att2 = [-0.03056 -0.207833 0.00728];
new_att2 = [-0.03056 -0.2398 0.009728];

old_length = sqrt(sum((att1 - old_att2).^2))
new_length = sqrt(sum((att1 - new_att2).^2))
