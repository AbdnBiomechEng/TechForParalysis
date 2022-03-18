function angles = estimate_shoulder_angles(hum_thory,hum_thorz,hum_thoryy)
% angles = [cy,cz,cx,sy,sz,sx,hy,hz,hyy] are the clavicular, scapular and
% humeral angles in the thorax, clavicle and scapula coordinate frames
% respectively (in radians)
%
% A regression model of the shoulder rhythm is used to estimate the
% clavicular and scapular angles from the humeral angles.
% The axial rotation of the clavicle is estimated separately.

% initial angles for clavicle and scapula, averaged among ?? subjects, from
% DSEM
C0=[-2.8573324e+001  1.1679804e+001  7.4820950e+000]*pi/180;
S0=[3.2047376e+001  1.0669100e+001 -1.0320259e+001]*pi/180;

% initial rotation matrices:
zero_clav = roty(C0(1))*rotz(C0(2))*rotx(C0(3));
zero_scap = roty(S0(1))*rotz(S0(2))*rotx(S0(3));

% linear regression model from de Groot and Brand 2001
[C_thory,C_thorz,S_thory,S_thorz,S_thorx] = shoulder_regression(C0,S0,hum_thory,hum_thorz);

% axial rotation of the clavicle: set to zero for now
C_thorx = 0;

% calculate rotation matrices
ThortoClav = roty(C_thory)*rotz(C_thorz)*rotx(C_thorx);
ThortoScap = roty(S_thory)*rotz(S_thorz)*rotx(S_thorx);
ThortoHum = roty(hum_thory)*rotz(hum_thorz)*roty(hum_thoryy);

% estimate axial rotation of the clavicle from initial matrices and
% scapular matrix
Ccorr=[zero_clav;ThortoClav];
Scorr=[zero_scap;ThortoScap];
[Cn,~]=axclav(Ccorr,Scorr);
ThortoClav = Cn(4:6,:);

% now find rotation matrices with respect to previous segment
ClavtoScap = ThortoClav'*ThortoScap;
ScaptoHum = ThortoScap'*ThortoHum;

% calculate Euler angles from these rotation matrices
Clavt = rotyzx(ThortoClav);
Scapc = rotyzx(ClavtoScap);
Hums = rotyzy(ScaptoHum);

angles = [Clavt Scapc Hums];

function [Cy,Cz,Sy,Sz,Sx] = shoulder_regression(C0,S0,Hy,Hz)
% estimates the clavicular and scapular angles from the humeral angles
% all angles with respect to thorax, in radians

% first change to degrees:
C0 = C0*180/pi;
S0 = S0*180/pi;
Hy = Hy*180/pi;
Hz = Hz*180/pi;

% regression values from de Groot and Brand 2001
% Clin Biomech Nov 16(9):735-43:
Cy = -4.983 + 0.120*Hy -0.242*Hz + 0.851*C0(1);
Cz = 3.917 - 0.046*Hy + 0.123*Hz + 0.493*C0(2);
Sy = -1.203 + 0.140*Hy -0.049*Hz + 0.901*S0(1);
Sz = 3.095 -0.079*Hy + 0.396*Hz + 0.414*S0(2);
Sx = 0.659 -0.028*Hy + 0.184*Hz +0.886*S0(3); 

% put angles back to radians:
Cy = Cy*pi/180;
Cz = Cz*pi/180;
Sy = Sy*pi/180;
Sz = Sz*pi/180;
Sx = Sx*pi/180;


 function angles = rotyzx(R)
% calculates the Euler angles around the y,z, and x axes
% from the rotation matrix R

z1 = asin(R(2,1));
sx =-R(2,3)/cos(z1);
cx = R(2,2)/cos(z1);
x1 = atan2(sx,cx);
sy =-R(3,1)/cos(z1);
cy = R(1,1)/cos(z1);
y1 = atan2(sy,cy);
if z1 >= 0,
  z2 = pi - z1;
else
  z2 = -pi - z1;
end
sx =-R(2,3)/cos(z2);
cx = R(2,2)/cos(z2);
x2 = atan2(sx,cx);
sy =-R(3,1)/cos(z2);
cy = R(1,1)/cos(z2);
y2 = atan2(sy,cy);
if (-pi/2 <= z1 && z1 <= pi/2)
  y=y1;
  z=z1;
  x=x1;
else
  y=y2;
  z=z2;
  x=x2;
end
angles = [y,z,x];


function angles=rotyzy(r)
% calculates the Euler angles around the y,z, and new y axes
% from the rotation matrix R

z1 = acos(r(2,2));
if (z1==0)
    y=acos(r(1,1));
	z=z1;
	ya=0.0;
	return;
end
sy = r(3,2)/sin(z1);
cy = -r(1,2)/sin(z1);
y1 = atan2(sy,cy);
sya = r(2,3)/sin(z1);
cya = r(2,1)/sin(z1);
ya1 = atan2(sya,cya);
z2 = -z1;
sy = r(3,2)/sin(z2);
cy = -r(1,2)/sin(z2);
y2 = atan2(sy,cy);
sya = r(2,3)/sin(z2);
cya = r(2,1)/sin(z2);
ya2 = atan2(sya,cya);
if (0 <= z1 && z1 <= pi)
   y = y1;
   z = z1;
   ya = ya1;
else
   y = y2;
   z = z2;
   ya = ya2;
end
% if z<(10*pi/180)
%     ya=y+ya;
%     y=0;
% end
angles = [y,z,ya];


function [Rx]=rotx(th)
% creates rotation matrix 
% for rotation of th radians around the x axis
Rx(1,1)=1;
Rx(2,2)=cos(th);
Rx(2,3)=-sin(th);
Rx(3,2)= sin(th);
Rx(3,3)= cos(th);


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


function[Cn,ROTr]=axclav(C,S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 02-08-91 %%
%                                                                       %
%  In de functie axclav wordt uit de positie matrices van clavicula en  %
%  scapula de axiale rotatie van de clavicula geschat en de nieuwe      %
%  positie matrix van de clavicula.                                     %
%                                                                       %
%       input : C (originele positiematrix van de clavicula)            %
%               S (positiematrix van de scapula)                        %
%       output: Cn (nieuwe positiematrix van de clavicula)              %
%               rot (rotaties van de clavicula rond de lokale assen)    %

%   --------------------------------------------------------------      %
%  The function axclav uses the original position matrices of the       %
%  clavicula and the the scapula to predict the axial rotation of the   %
%   clavicula and the new position matrix of the clavicula.             %
%
%       input : C (original position matrix of the clavicula)           %
%               S (position matrix of the scapula)                      %
%       output : Cn (new position matrix of the clavicula)              %
%                rot (rotations of the clavicula about the local axes)  %
%                                                                       %
%  Please note that the first frame of data C(1:3,1:3) and S(1:3,1:3)   %
%  is at the resting (zero) position                                    %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ruststand van clavicula en scapula:
% Rest position of the clavicula and scapula:
  Co = C(1:3,1:3); So = S(1:3,1:3);

% Relatie clavicula-scapula in ruststand:
% clavicula-scapula relationship in rest:
  Mp = Co'*So;

% Bepalen van de rotaties van de clavicula rond de lokale assen t.o.v.
% thorax LCS:
% Determine the rotations of the clavicula about the local axis w.r.t.
% thorax LCS:
[n,m]=size(C);nDATA = n/3;

r=[];
for i=1:nDATA
  Rc=Co'*C(3*i-2:3*i,1:3);
  r = [r;Rc];
end  

% r = rot2a(C,0);

% Bepalen van de eulerrotaties (y,z,x):
% Determine the eulerrotations (y,z,x):
  rot1=rotyzx(r(1:3,1:3));
  rot2=rotyzx(r(4:6,1:3));
  ROTr = [rot1;rot2];
  [rij,kolom]=size(ROTr);

  BETA = ROTr(:,1); GAMMA = ROTr(:,2); ALPHA = ROTr(:,3);

for I=1:rij
  beta=BETA(I);gamma=GAMMA(I);          % Benodigde hoeken
                                        % Required angles
  RS  = S(I*3-2:I*3,1:3);                 % Scapula stand i
                                          % position/orientation of the
                                          % scapula i
%____________________________________________________________________

% (Berekening van axiale rotatie volgens axrot (vd Helm))
% Calculation of the axial rotations according to axrot (vd Helm)
 
    Ma = Co*roty(beta)*rotz(gamma);
    alpha=0;
    Sposd = Ma*rotx(alpha)*Mp;
    Spos=RS;
    Emat = Sposd'*Spos;
    E = acos(Emat(1,1)) + acos(Emat(2,2)) + acos(Emat(3,3));
    SSQ = E*E;
    SSQo = 0;
    while abs(SSQ-SSQo)>0.001,
      dalpha = max(0.1,alpha)*sqrt(eps);
      Sposdd = Ma*rotx(alpha + dalpha)*Mp;
      Ematd = Sposdd'*Spos;
      Ed = acos(Ematd(1,1)) + acos(Ematd(2,2)) + acos(Ematd(3,3));
      dEdalpha = (Ed - E)/dalpha;
      V = E/dEdalpha;
      d=2;
      alpha0=alpha;
      SSQo = SSQ;
      while (SSQo <= SSQ) && ((abs(SSQo - SSQ) > 0.001)||(d==2)),
	d = d/2;
	alpha1 = alpha0 - d*V;
	Sposd = Ma*rotx(alpha1)*Mp;
	Emat = Sposd'*Spos;
	E = acos(Emat(1,1)) + acos(Emat(2,2)) + acos(Emat(3,3));
	SSQ=E*E;
      end
      alpha=alpha1;

    end

    ALPHA(I) = alpha; 
end

  ROTr(1:rij,3)=ALPHA;
  ROTg(1:rij,3)=ALPHA*180/pi;

% Berekenen van de nieuwe positiematrix van de clavicula :
% (volgens Ci = Co * R waarbij R = (roty*rotz*rotx)
% Calculation of the new position matrix of the clavicula:
% according to Ci = Co * R with R = (roty*rotz*rotx)

for I=1:rij
  Ry = roty(BETA(I));
  Rz = rotz(GAMMA(I));
  Rx = rotx(ALPHA(I));
  Cn(I*3-2:I*3,1:3) = Co*Ry*Rz*Rx;

end

