all_thor = zeros(10,1);
all_scap = zeros(10,2);
for i=1:20
    Result = das3_optimize_hang('random',['random_init_' num2str(i)]);
    all_thor(i) = Result.FGHcontact;
    all_scap(i,:) = Result.Fscap;
end