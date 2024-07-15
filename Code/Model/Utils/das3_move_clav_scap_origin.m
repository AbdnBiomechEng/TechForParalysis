%% Read in model

% import OpenSim namespace
import org.opensim.modeling.*;

% Create the Original OpenSim model from a .osim file
Model1 = Model('das3.osim');
Model1.initSystem;

% Create a copy of the original OpenSim model for the Modified Model
Model2 = Model(Model1);
Model2.initSystem;

% Rename the modified Model so that it comes up with a different name in
% the GUI navigator
Model2.setName('DAS3_clav_scap_orig');

%% Get conversion factors for new location of SC and AC body origins

% Get the set of joints
JointSet = Model1.getJointSet();

% Find the SC joint and get position of joint in body origin
SCjoint = JointSet.get('sc1');

child_frame=SCjoint.get_frames(1);
p = child_frame.get_translation();
SC_in_SCorigin = [p.get(0) p.get(1) p.get(2)];

% Find the AC joint and get position of joint in body origin
ACjoint = JointSet.get('ac1');

child_frame=ACjoint.get_frames(1);
p = child_frame.get_translation();
AC_in_ACorigin = [p.get(0) p.get(1) p.get(2)];


%% Joint positions

% In new model, CS and AC are the origins of the clavicle and scapular
% frames, so the child frame translation is [0,0,0]
Joint2 = Model2.getJointSet();

SCjoint2 = Joint2.get('sc1');
SCjoint2.get_frames(1).set_translation(Vec3(0,0,0));
SCjoint2 = Joint2.get('sc2');
SCjoint2.get_frames(0).set_translation(Vec3(0,0,0));
SCjoint2.get_frames(1).set_translation(Vec3(0,0,0));
SCjoint2 = Joint2.get('sc3');
SCjoint2.get_frames(0).set_translation(Vec3(0,0,0));
SCjoint2.get_frames(1).set_translation(Vec3(0,0,0));

ACjoint2 = Joint2.get('ac1');
p = ACjoint2.get_frames(0).get_translation();
AC_origin = [p.get(0) p.get(1) p.get(2)];
p = AC_origin - SC_in_SCorigin;
ACjoint2.get_frames(0).set_translation(Vec3(p(1),p(2),p(3)));
ACjoint2.get_frames(1).set_translation(Vec3(0,0,0));
ACjoint2 = Joint2.get('ac2');
ACjoint2.get_frames(0).set_translation(Vec3(0,0,0));
ACjoint2.get_frames(1).set_translation(Vec3(0,0,0));
ACjoint2 = Joint2.get('ac3');
ACjoint2.get_frames(0).set_translation(Vec3(0,0,0));
ACjoint2.get_frames(1).set_translation(Vec3(0,0,0));

GHjoint2 = Joint2.get('gh1');
p = GHjoint2.get_frames(0).get_translation();
GH_origin = [p.get(0) p.get(1) p.get(2)];
p = GH_origin - AC_in_ACorigin;
GHjoint2.get_frames(0).set_translation(Vec3(p(1),p(2),p(3)));

%% Mass Centres

% We also need to update the mass centres of the clavicle and scapule
Body2 = Model2.getBodySet();

clavicle = Body2.get('clavicle_r');
p = clavicle.getMassCenter;
clav_mass_centre = [p.get(0) p.get(1) p.get(2)];
p = clav_mass_centre - SC_in_SCorigin;
clavicle.set_mass_center(Vec3(p(1),p(2),p(3)));

scapula = Body2.get('scapula_r');
p = scapula.getMassCenter;
scap_mass_centre = [p.get(0) p.get(1) p.get(2)];
p = scap_mass_centre - AC_in_ACorigin;
scapula.set_mass_center(Vec3(p(1),p(2),p(3)));

%% Muscle attachment points

% Get the set of muscles that are in the original model
MuscleSet = Model1.getMuscles(); 
nMuscles = MuscleSet.getSize();

% Get the set of muscles that are in the new model
% (Should be the same as the original at this point.)
Muscles2 = Model2.getMuscles();

% loop through muscles and change attachment points as needed
for imus = 0:nMuscles-1
        
    currentMuscle = MuscleSet.get(imus);
    newMuscle = Muscles2.get(imus);

    % get muscle path
    muspath = currentMuscle.getGeometryPath();
    PtSet = muspath.getPathPointSet();
    nPts = PtSet.getSize();

    newmuspath = newMuscle.getGeometryPath();
    newPtSet = newmuspath.getPathPointSet();

    for ipt = 0:nPts-1
        if strcmp(PtSet.get(ipt).getBodyName(),'clavicle_r')
            point = PathPoint.safeDownCast(PtSet.get(ipt));
            p = point.get_location();
            pt_in_prev_frame = [p.get(0) p.get(1) p.get(2)];
            p = pt_in_prev_frame - SC_in_SCorigin;
            newpoint = PathPoint.safeDownCast(newPtSet.get(ipt));
            newpoint.set_location(Vec3(p(1),p(2),p(3)));

        elseif strcmp(PtSet.get(ipt).getBodyName(),'scapula_r')
            point = PathPoint.safeDownCast(PtSet.get(ipt));
            p = point.get_location();
            pt_in_prev_frame = [p.get(0) p.get(1) p.get(2)];
            p = pt_in_prev_frame - AC_in_ACorigin;
            newpoint = PathPoint.safeDownCast(newPtSet.get(ipt));
            newpoint.set_location(Vec3(p(1),p(2),p(3)));
        end
    end

end

%% Markers

% Get the set of markers that are in the original model
MarkerSet = Model1.getMarkerSet(); 
nMarkers = MarkerSet.getSize();

% Get the set of markers that are in the new model
% (Should be the same as the original at this point.)
Markers2 = Model2.getMarkerSet();

% loop through markers and change location as needed
for imark = 0:nMarkers-1
        
    currentMarker = MarkerSet.get(imark);
    newMarker = Markers2.get(imark);

    if strcmp(currentMarker.getParentFrame(),'clavicle_r')
        p = currentMarker.get_location();
        pt_in_prev_frame = [p.get(0) p.get(1) p.get(2)];
        p = pt_in_prev_frame - SC_in_SCorigin;
        newMarker.set_location(Vec3(p(1),p(2),p(3)));

    elseif strcmp(currentMarker.getParentFrame(),'scapula_r')
        p = currentMarker.get_location();
        pt_in_prev_frame = [p.get(0) p.get(1) p.get(2)];
        p = pt_in_prev_frame - AC_in_ACorigin;
        newMarker.set_location(Vec3(p(1),p(2),p(3)));
    end

end


%% Save the updated model to an OSIM xml file
Model2.print('das3_clav_scap_orig.osim');
disp('The new model has been saved at das3_clav_scap_orig.osim');


