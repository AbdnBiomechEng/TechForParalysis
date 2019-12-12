% We want to read trc file and starting
%from line (row) 7, then
%extract data from coloumn 3, repeated three times.
% look out for the names C7, PX, T8, IJ, SC, AC, AI, TS, AA, GH, LE, ME,
% US, RS. then at the end the code should output trc file with the
% transformed data.

K=readtrc('NewSession101.trc');
%K=OriginalMatrix';   %loading the data
[m,n]=size(OriginalMatrix);
%n12 = length (K);

%Rotation parameter & transformation
% x to z, y to x, z to y 
rx1 = pi/2;       %Rotation value (rad) in x - change as needed
ry1 = pi;       %Rotation value (rad) in y - change as needed
rz1 = pi/2;       %Rotation value (rad) in z - change as needed

Rx = [1 0 0;
    0 cos(rx1) -sin(rx1);
    0 sin(rx1) cos(rx1)];

Ry = [cos(ry1) 0 sin(ry1);
      0 1 0;
      -sin(ry1) 0 cos(ry1)];
  
Rz = [cos(rz1) -sin(rz1) 0;
      sin(rz1) cos(rz1) 0;
      0 0 1];

R = Rx*Ry*Rz;                    %Rotation matrix 

%Transformed matrix
n3=n/3; %no. of markers
Rrepeated=repmat(R,n/3);
M1 = Rrepeated*K;   %Apply rotation parameters to transform
M1=M1';

%Translation parameter & transformation
%539.70264	-14.28845	1435.79712 (x,y,z - original)
Tx1 = M1(1,1);     %Translation value in x - change as needed
Ty1 = M1(1,2);     %Translation value in y - change as needed
Tz1 = M1(1,3);     %Translation value in z - change as needed
T1 = [Tx1; Ty1; Tz1];
%Trepeated=repmat(T1,n/3,n/3);
%TransformedData=M1-Trepeated; %Apply rotation and translation


