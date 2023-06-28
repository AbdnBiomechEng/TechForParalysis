% Script to add ground segment to model struct
% The ground segment is not longer included in the new opensim format, but
% the realtime model requires it

modelparams = load('simplified_model_struct.mat'); 
model = modelparams.model;
model.nSegments = 14;
segments = model.segments;

ground = load('ground_segment.mat'); 
model.segments{1} = ground.ground_segment;
for i=2:14
    model.segments{i} = segments{i-1};
end
save('simplified_model_struct','model');