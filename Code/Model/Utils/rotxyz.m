function [x,y,z] = rotxyz(R)
% calculates the Euler angles around the x,y, and z axes
% from the rotation matrix R

y1 = asin(R(1,3));
sz = -R(1,2)/cos(y1);
cz =  R(1,1)/cos(y1);
z1 = atan2(sz,cz);
sx = -R(2,3)/cos(y1);
cx =  R(3,3)/cos(y1);
x1 = atan2(sx,cx);
if y1>=0 
  y2 = pi - y1;
else
  y2 = -pi -y1;
end
sz = -R(1,2)/cos(y2);
cz =  R(1,1)/cos(y2);
z2 = atan2(sz,cz);
sx = -R(2,3)/cos(y2);
cx =  R(3,3)/cos(y2);
x2 = atan2(sx,cx);
if (-pi/2 <= y1 & y1 <= pi/2)
  y=y1;
  z=z1;
  x=x1;
else
  y=y2;
  z=z2;
  x=x2;
end