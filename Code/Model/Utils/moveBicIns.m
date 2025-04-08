% moveBicIns.m
% Script to estimate updated insertion positions of biceps - moving from
% one insertion point to one proximal on the tuberosity to one distal
% according to van Bekerom 2014
% https://link.springer.com/article/10.1007/s00167-014-3322-9


cyl = [-1.6165 0.18302 1.4708]; % orientation of wrapping cylinder in radius frame (Euler angles)

Rot_cyl2rad = eul2rotm(cyl, "XYZ"); % rotation matrix to convert cylinder frame to radius frame
Bicl_transl_cyl = [0 0 0.00513]'; % translation of BIC_L in cyl frame
Bicl_transl_rad = Rot_cyl2rad*Bicl_transl_cyl; % translation in radius frame

Bicb_transl_cyl = [0 0 -0.00513]'; % tranlation of BIC_B in cyl frame
Bicb_transl_rad = Rot_cyl2rad*Bicb_transl_cyl; % translation in radius frame

bic_ins = [-0.00839 -0.036287 -0.0061803]'; % old single insertion point

bicl_ins = (bic_ins + Bicl_transl_rad)' % new BicL insertion
bicb_ins = (bic_ins + Bicb_transl_rad)' % new BicB insertion

