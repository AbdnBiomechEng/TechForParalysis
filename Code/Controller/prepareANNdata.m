% Prepare ANN training data
high = load('ANNtraining/high_reach');
side = load('ANNtraining/side_reach');
across = load('ANNtraining/reach_across');
mouth = load('ANNtraining/reach_to_mouth');
mouth2 = load('ANNtraining/reach_to_mouth2');
side2 = load('ANNtraining/side_reach2');
side3 = load('ANNtraining/side_reach3');

nmus = high.Result.nmus;

musclenames = cell(nmus,1);
for imus=1:nmus
    musclenames{imus} = high.Result.muscles{imus}.name;
end

include_mus = [12 13 14 21];
tric = 29;
%include_mus = [12 13 14 21 22];

% inputs: bic and delt
inputs_orig = [
side2.Result.u(include_mus,:),...
mouth.Result.u(include_mus,:),...
high.Result.u(include_mus,:),...
side.Result.u(include_mus,:),...
mouth2.Result.u(include_mus,:),...
across.Result.u(include_mus,:)];

% output: tric
output_orig = [
    side2.Result.u(tric,:),...
    mouth.Result.u(tric,:),...
    high.Result.u(tric,:),...
    side.Result.u(tric,:),...
    mouth2.Result.u(tric,:),...
    across.Result.u(tric,:)];

timevec_orig = 3:3:(6*120);
timevec = 0.3:0.3:(6*120);

for i = 1:length(include_mus)
    inputs_orig(i,:) = smooth(inputs_orig(i,:));
    inputs(i,:) = interp1(timevec_orig,inputs_orig(i,:),timevec,'spline');

end
output_orig = smooth(output_orig)';
output = interp1(timevec_orig,output_orig,timevec,'spline');

timevec_orig = 3:3:120;
timevec = 0.3:0.3:120;

inputs_test_orig = side3.Result.u(include_mus,:);
output_test_orig = side3.Result.u(tric,:);

for i = 1:length(include_mus)
    inputs_test_orig(i,:) = smooth(inputs_test_orig(i,:));
    inputs_test(i,:) = interp1(timevec_orig,inputs_test_orig(i,:),timevec,'spline');

end
output_test_orig = smooth(output_test_orig)';
output_test = interp1(timevec_orig,output_test_orig,timevec,'spline');

