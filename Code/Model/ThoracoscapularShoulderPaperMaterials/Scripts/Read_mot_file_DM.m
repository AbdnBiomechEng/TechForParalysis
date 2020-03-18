motData = dlmread('S08_IK_WC3.mot', '\t', 11, 0);
import org.opensim.modeling.*
trialsForCMC = dir('*.mot'); %reads all current .mot files in your working directory
nTrials = size(trialsForCMC);
trial = 1;
motFile = trialsForCMC(trial).name;
coordinateSto = Storage(fullfile(motFile));
%coordinates = model.getCoordinateSet();
Time = ArrayDouble();
nFrames = coordinateSto.getTimeColumn(Time);
nCoordinates = coordinateSto.getColumnLabels().getSize()-1;
for dof = 0:nCoordinates-1
    coordvalue = ArrayDouble();
    for fr = 0:nFrames-1
        coordinateSto.getDataColumn(dof, coordvalue);
        motData(fr+1, dof+1) = coordvalue.getitem(fr);
    end
end