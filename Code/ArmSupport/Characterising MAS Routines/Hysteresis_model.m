%% Add hysteresis to the model of the relationship 
%  between Spring Setting, Arm Height and Force

% The bulk of this routine is identical to the first sections of the
% Force_Height_Spring_model routine, but the hysteresis effect is added at the
% end
%%
clear
close all
%%  Define the additional variables needed to include hysteresis

% hyst_range_asc and hyst_range_desc:  These set the change in height needed to
% completely reverse the frictional effect when the arm changes direction of
% travel
% exFfric: The 'system friction coefficient', determined in
% Force_Height_Spring_model, used to calculate the frictional force to be added
% to or subtracted from the average force, to model the hysteresis effect

hyst_range_asc=40;  % mm  when arm starts moving up
hyst_range_desc=65; % mm  when arm starts moving down
exFfric=0.21;       

%%  Import the data

fd_orig=load('fd_orig.mat');
data=fd_orig.data;
numdatasets=size(data,2);
fsamp=100;
tstep=1/fsamp;

%%  Reset linestyle order to default 

set(groot,'defaultAxesColorOrder', 'remove');
set(groot,'defaultAxesLineStyleOrder', 'remove');
%  Set linestyles to my preferred colourblind styles..
prefcolours=[1 0 0;0 0 1;0 1 0;0 0 0];
set(groot,'defaultAxesColorOrder',prefcolours,...
     'defaultAxesLineStyleOrder','-|--|:')

%% Check how starting height and force zero offset vary with spring setting

for n=1:numdatasets
        datacell=data(1,n);
        dataset=datacell{1,1};
        startheight(n)=dataset(1,1);
        zerooffset(n)=dataset(1,2);
        springsetting(n)=data{3,n};
  end

 figure;
 subplot(2,1,1)
    scatter(springsetting, startheight)  
    title('Starting Height vs. Spring Setting');
    ylabel('Starting Height, mm');
    %xlabel('Spring Setting');
 subplot(2,1,2)
    scatter(springsetting, zerooffset)  
    title('Force Zero Offset vs. Spring Setting');
    ylabel('Force Zero Offset, N');
    xlabel('Spring Setting');
    
    %% Adjust zero offset for outliers
 
% Assume incorrect zero offset for two of the Spring setting 6 tests 
% and the spring setting 15 test: 
% Assume force already being applied when test began
% Change these zero offsets to the average value from the other tests

   avgoffset(1:14)=zerooffset(1:14);
   avgoffset(15:19)=zerooffset(17:21);
   avgoffset(20)=zerooffset(23);
   avgoffsetfinal=mean(avgoffset)
   data{1,15}(1,2)=avgoffsetfinal;
   data{1,16}(1,2)=avgoffsetfinal;
   data{1,22}(1,2)=avgoffsetfinal; 
   
   %% Make Initial Plots
   
   % Initial plots of Force vs. Displacement and 
   % Force and Displacement vs. Time
   
   %reset linestyle order to default 
    set(groot,'defaultAxesColorOrder', 'remove');
    set(groot,'defaultAxesLineStyleOrder', 'remove');
    %set linestyles to my preferred colourblind styles..
    set(groot,'defaultAxesColorOrder',prefcolours,...
     'defaultAxesLineStyleOrder','-|--|:')
 
 % Set cut off for low pass filter
 
    f_cutoff = .5;
 
% Loop to process data
    
for n=1:numdatasets
    %Collect set of data
    datacell=data(1,n);
    testsetting=data(2,n);
    dataset=datacell{1,1};
    %Create time vector
    t=linspace(1,size(dataset,1),size(dataset,1));
    t=t*tstep;
    %Correct the force values with the zero offset
    foffset=dataset(1,2);
    dataset(:,2)=dataset(:,2)-foffset;
    %Adjust the height values to make them displacement values relative to
    %the starting height
    yoffset=dataset(1,1);
    dataset(:,1)=dataset(:,1)-yoffset; 
    
    % Build and apply low-pass filter

    % Build low-pass Butterworth filter
    [b,a] = butter(2,f_cutoff/(fsamp/2),'low');  
    % Apply filter to data                                                                                         
    force_filt = filtfilt(b,a,dataset(:,2));
     height_filt = filtfilt(b,a,dataset(:,1));
    
    figure;
    subplot(2,1,1);
    plot(dataset(:,1),dataset(:,2))
    hold on
    plot(height_filt,force_filt,'-k')
    title([testsetting ' Force vs. Displacement']);
    ylabel('Force, N)');
    xlabel('Displacement, mm');
    subplot(2,1,2);
    title(['Variation with time']);
    xlabel('time, secs');
    yyaxis right
    plot(t,dataset(:,1))
    hold on
    plot(t,height_filt,':b')
    hold on
    ylabel('Displacement, mm');
    yyaxis left
    plot(t,dataset(:,2))
    hold on
    plot(t,force_filt,'-k')
    ylabel('Force, N');
    hold off
end
%%  Remove two outlier points from dataset 5
   idx=find(~all(data{1,5},2));
   data{1,5}(idx,:)=[];
%  Checked by running previous section again: Outliers have been removed
%% Make plots for each spring setting, including hysteresis

% Some plots to compare results at each spring setting
% These also include the predicted force values using Equation 1 to generate the
% average force, based on the height, and then the linear proportional model of
% hysteresis used in the final forward dynamic simulations


%reset linestyle order to default
set(groot,'defaultAxesColorOrder', 'remove');
set(groot,'defaultAxesLineStyleOrder', 'remove');
%set linestyles to my preferred colourblind styles..
set(groot,'defaultAxesColorOrder',prefcolours,...
     'defaultAxesLineStyleOrder','-|--|:')
 
f_cutoff = .5;

% Set the working height range
minheight=500;
maxheight=1100;

%  Loop which plots data for each spring setting

for s=1:15
% just plot spring setting 4 on oscillations for report
% for s=4:4   
    legendsetting=[];
    figure;
   
   for n=1:numdatasets
   % just plot oscillations test for report
   % for n=11:11
    springsetting=data{3,n};
    if springsetting==s
        datacell=data(1,n);
        testsetting=data(2,n);
        legendsetting=[legendsetting,testsetting,'Modelled'];
        dataset=datacell{1,1};
        t=linspace(1,size(dataset,1),size(dataset,1));
        t=t*tstep;
        foffset=dataset(1,2);
        dataset(:,2)=dataset(:,2)-foffset;
        
        % Build low-pass Butterworth filter
        [b,a] = butter(2,f_cutoff/(fsamp/2),'low');  
        % Apply filter to data                                                                                         
        force_filt = filtfilt(b,a,dataset(:,2));
        height_filt = filtfilt(b,a,dataset(:,1));
        
         % filter out values where height is ouside working range
        idxlow=find(height_filt(:,1)<minheight);
        force_filt(idxlow,:)=[0];
        height_filt(idxlow,:)=[0];
        idxhigh=find(height_filt(:,1)>maxheight);
        force_filt(idxhigh,:)=[0];
        height_filt(idxhigh,:)=[0];
        force_filt=force_filt(all(force_filt,2),:);
        height_filt=height_filt(all(height_filt,2),:);  
        
        hold on
        % Only the low-pass filtered data is plotted
        plot(height_filt,force_filt)
        
        % create the modelled average force, and adjust for hysteresis
        % Force = -5.778*SpringSetting +0.0073*(Height-500) -3.854
        
        force_model=zeros(1,length(height_filt));
        force_avg=zeros(1,length(height_filt));
        for h=1:length(height_filt)
              force_avg(h)=-5.778*s+0.0073*(height_filt(h)-500)-3.854;
           if h>1
               force_avg_change=force_avg(h)-force_avg(h-1);
               exFzchange=height_filt(h)-height_filt(h-1);
              if exFzchange<0  % arm descending
                force_model(h)=max(force_model(h-1)+force_avg_change...
                    -(exFzchange/hyst_range_desc)*force_avg(h)...
                    *exFfric,force_avg(h)*(1+exFfric));
            else if exFzchange>0 % arm ascending
                force_model(h)=min(force_model(h-1)+force_avg_change...
                    -(exFzchange/hyst_range_asc)*force_avg(h)...
                    *exFfric,force_avg(h)*(1-exFfric));
                else
                    force_model(h)=force_avg(h);
                end
              end
           else
             force_model(h)=force_avg(h);
           end            
        end
        plot(height_filt,force_model)
        
    else
    end
  end
            title(['Spring setting ',int2str(s),' Force vs. Height']);
            ylabel('Force, N');
            xlabel('Height, mm');
            legend(legendsetting);
    hold off
end

