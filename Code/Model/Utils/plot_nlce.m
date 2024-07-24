% script to plot normalised fibre lengths saved using das3_check_lceopt.m

load('all_lengths_lce.mat');

% choose which muscle elements to plot
allmus = 89;

% it addes a horizontal line at 1 (where fibre length = optimal fibre
% length)
figure; 
for imus=allmus
    subplot(length(allmus),1,imus-allmus(1)+1); plot(n_lce(imus,:)); hold on; yline(1); title(musnames{imus},'Interpreter', 'none');
end