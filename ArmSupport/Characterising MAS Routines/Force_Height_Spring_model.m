%% Create a model of the relationship between 
%  Spring Setting, Arm Height and Force
%%  Import the data

clearclose all
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
     disp_filt = filtfilt(b,a,dataset(:,1));
    
    figure;
    subplot(2,1,1);
    plot(dataset(:,1),dataset(:,2))
    hold on
    plot(disp_filt,force_filt,'-k')
    title([testsetting ' Force vs. Displacement']);
    ylabel('Force, N)');
    xlabel('Displacement, mm');
    subplot(2,1,2);
    title(['Variation with time']);
    xlabel('time, secs');
    yyaxis right
    plot(t,dataset(:,1))
    hold on
    plot(t,disp_filt,':b')
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
%% Make plots for each spring setting

% Some plots to compare results at each spring setting

%reset linestyle order to default
set(groot,'defaultAxesColorOrder', 'remove');
set(groot,'defaultAxesLineStyleOrder', 'remove');
%set linestyles to my preferred colourblind styles..
set(groot,'defaultAxesColorOrder',prefcolours,...
     'defaultAxesLineStyleOrder','-|--|:')
 
f_cutoff = .5;
 
%  Loop which plots data for each spring setting

for s=1:15
    legendsetting=[];
    figure;
    for n=1:numdatasets
    springsetting=data{3,n};
    if springsetting==s
        datacell=data(1,n);
        testsetting=data(2,n);
        legendsetting=[legendsetting,testsetting];
        dataset=datacell{1,1};
        t=linspace(1,size(dataset,1),size(dataset,1));
        t=t*tstep;
        foffset=dataset(1,2);
        dataset(:,2)=dataset(:,2)-foffset;
        yoffset=dataset(1,1);
        dataset(:,1)=dataset(:,1)-yoffset;
        
        % Build low-pass Butterworth filter
        [b,a] = butter(2,f_cutoff/(fsamp/2),'low');  
        % Apply filter to data                                                                                         
        force_filt = filtfilt(b,a,dataset(:,2));
        disp_filt = filtfilt(b,a,dataset(:,1));
        
        %hold on
        % Original data not plotted as congests the figures
        % plot(dataset(:,1),dataset(:,2));
        hold on
        % Only the low-pass fitered data is plotted
        plot(disp_filt,force_filt)
    else
    end
  end
            title(['Spring setting ',int2str(s),' Force vs. Displacement']);
            ylabel('Force, N)');
            xlabel('Displacement, mm');
            legend(legendsetting);
    hold off
end

%%  Note on processing data

   % The hysteresis evident on the plots will complicate adding the force
   % to the model.  Therefore the first step will be to ignore the
   % hysteresis and use a single value to represent the force at a given
   % height for each spring setting.
   
   % A first hypothesis is that the hysteresis is due to friction in the
   % pulley system, and if this could be removed, the single representative
   % value will be the average of the raising and lowering values for a
   % given height and spring setting.
   
   % There also appears to be some 'slack' in the system at small
   % displacements once the spring setting is above 2.  There are also some
   % small variations at the lower end of travel, for example if the
   % applied force was increased after the mechanism reached its lower
   % stop.
   
   % Therefore, the data will be analysed over a height range of 500 -
   % 1100mm, giving a practical working window of 600mm vertical change for
   % the user.
   
   % Also, the hysteresis loops created by the oscillating tests affect the
   % fitting of linear regression lines to the ascending and descending
   % data at a later stage.  Therefore these datasets are excluded at this
   % stage
   %%  Retain one down-up dataset for each spring setting
   keep = [5 8 10 12 14 17 18 19 21 22];
   data=data(:,keep);
   numdatasets=size(data,2);
      %%  Extract and plot filtered force data within chosen height range
   
   % Repeat previous step, but plot height rather than displacement, and
   % only include values with heights in the range.
   %  Also create a cell array with all the values for a given spring
   %  setting in the height range in each cell
  
%reset linestyle order to default 
set(groot,'defaultAxesColorOrder', 'remove');
set(groot,'defaultAxesLineStyleOrder', 'remove');
%set linestyles to my preferred colourblind styles..
set(groot,'defaultAxesColorOrder',prefcolours,...
     'defaultAxesLineStyleOrder','-|--|:')
 
f_cutoff = .5;
% Set the height range
minheight=500;
maxheight=1100;
% Cell array to contain data within the height range
rangeout=cell(1,15);

%  Loop to extract data at each spring setting as before

for s=1:15
    
    for n=1:numdatasets
    springsetting=data{3,n};
    
    if springsetting==s
        datacell=data(1,n);
        testsetting=data(2,n);
        dataset=datacell{1,1};
        foffset=dataset(1,2);
        dataset(:,2)=dataset(:,2)-foffset;
        
       % Build low-pass Butterworth filter
        [b,a] = butter(2,f_cutoff/(fsamp/2),'low');  
        % Apply filter to data                                                                                         
        force_filt = filtfilt(b,a,dataset(:,2));
        height_filt = filtfilt(b,a,dataset(:,1));
        
        % filter out values where height is ouside range
        idxlow=find(height_filt(:,1)<minheight);
        force_filt(idxlow,:)=[0];
        height_filt(idxlow,:)=[0];
        idxhigh=find(height_filt(:,1)>maxheight);
        force_filt(idxhigh,:)=[0];
        height_filt(idxhigh,:)=[0];
        force_filt=force_filt(all(force_filt,2),:);
        height_filt=height_filt(all(height_filt,2),:);  

        rangeout{1,s}=[rangeout{1,s}; height_filt force_filt];
       
    else
    end
         
    end
    
    if ~isempty(rangeout{1,s})
        figure;
        out1=rangeout{1,s}(:,1);
        out2=rangeout{1,s}(:,2);
        plot(out1,out2)
        titleset=['Force vs. Height over chosen range, Spring Setting ',int2str(s)];
        title(titleset);
        ylabel('Force, N');
        xlabel('Height, mm');
   
    else
    end
    
end
%%  Linear regression analysis of chosen data

% First order LS regression lines added.  Adding a single line to each data
% set is found to give odd results in some cases.  So the data is split
% into two parts, representing the ascending and descending movement of the
% MAS arm, and linear regression lines added for each part.  Then a new
% central regression line is found by averaging the coefficients for the
% ascending and descending parts.  This generates a representative line in
% all cases.

%reset linestyle order to default 
set(groot,'defaultAxesColorOrder', 'remove');
set(groot,'defaultAxesLineStyleOrder', 'remove');
%set linestyles to my preferred colourblind styles..
set(groot,'defaultAxesColorOrder',prefcolours,...
     'defaultAxesLineStyleOrder','-|:|--')

regressout=zeros(15,2);  % Array to collect regression coefficients
rsqout=[];               % Array to collect R values 
fresidout=[];            % Array to collect residuals
ascending=cell(1,15);    % Array to collect ascending values
descending=cell(1,15);   % Array to collect descending values
hysteresis=[];           % Array to collect some data on hysteresis

for cc=1:15
    if ~isempty(rangeout{1,cc})
        datacell=rangeout{1,cc};
        height=datacell(:,1);
        force=datacell(:,2);
        % fit linear regression line to data
        p=polyfit(height,force,1);
        % find predicted value
        calcforce=polyval(p,height);
        
        % Find R value to gauge fit
        % find residual
        fresid=force-calcforce;
        % find R value
        SSresid=sum(fresid.^2);
        SStotal = (length(force)-1) * var(force);
        rsq = 1 - SSresid/SStotal;
        %plot results        
        figure;
        %subplot(2,1,1)
        plot(height,force)
        hold on
        plot (height,calcforce)
        
        % use an alternative method for comparison
        % split data into ascending and descending parts
        idxasc=find(fresid>0);
        idxdesc=find(fresid<=0);
        heightasc=height(idxasc);
        forceasc=force(idxasc);
        heightdesc=height(idxdesc);
        forcedesc=force(idxdesc);
        ascending{1,cc}=[heightasc forceasc];
        descending{1,cc}=[heightdesc forcedesc];
        % create separate regression lines for ascending and descending
        pasc=polyfit(heightasc,forceasc,1);
        calcforceasc=polyval(pasc,heightasc);
        pdesc=polyfit(heightdesc,forcedesc,1);
        calcforcedesc=polyval(pdesc,heightdesc);
        %plot those lines
        hold on
        plot(heightasc,calcforceasc)
        hold on
        plot(heightdesc,calcforcedesc)
        
        % Find R values for ascending and descending regression line
        % Ascending:
        % find residual
        fresidasc=forceasc-calcforceasc;
        % find R value
        SSresidasc=sum(fresidasc.^2);
        SStotalasc = (length(forceasc)-1) * var(forceasc);
        rsqasc = 1 - SSresidasc/SStotalasc;
        % Descending:
        % find residual
        fresiddesc=forcedesc-calcforcedesc;
        % find R value
        SSresiddesc=sum(fresiddesc.^2);
        SStotaldesc = (length(forcedesc)-1) * var(forcedesc);
        rsqdesc = 1 - SSresiddesc/SStotaldesc;
        
        %  create a new set of coefficients for a central line by combining
        %  the ascending and descending coefficients
        pcentral=(pasc+pdesc)/2;
        calcforcecentral=polyval(pcentral,height);
        hold on
        plot(height,calcforcecentral)
        ylabel('Force, N');
        xlabel('Height, mm');
        legend('Original Data','LS All Data','LS Ascending','LS Descending','LS Central')
        title(['Spring setting ',int2str(cc),' Force vs. Displacement']);
        
        % save coefficients
        regressout(cc,:)=pcentral';
        
        % Save some info on hysteresis at max and min height, to see how it
        % varies with spring setting
        
        hysteresis=[hysteresis;cc min(calcforceasc) min(calcforcecentral) min(calcforcedesc) max(calcforceasc) max(calcforcecentral) max(calcforcedesc)]; 
        
        %  Find R value for central line (out of interest)
        % find residual
        fresidcent=force-calcforcecentral;
        % find R value
        SSresidcent=sum(fresidcent.^2);
        rsqcent = 1 - SSresidcent/SStotal;
        rsqout=[rsqout;cc rsq rsqcent rsqasc rsqdesc];

        
        hold off
        %  This plot not used for final report
        %subplot(2,1,2)
        %scatter(height,fresid,1,'.')
        % hold on
        %scatter(height,fresidcent,1,'.')
        %scatter(heightasc,fresidasc,1,'.')
        %scatter(heightdesc,fresiddesc,1,'.')
        %legend('Residual All data','Residuals Central','Residuals Ascending','Residuals Descending')
        %xlabel('Height');
        %ylabel('Residual Value');
        %title(['Spring setting ',int2str(cc),' Residuals']);
    else 
    end
end

figure;
scatter(rsqout(:,1),rsqout(:,2))
hold on
scatter(rsqout(:,1),rsqout(:,3))
scatter(rsqout(:,1),rsqout(:,4))
scatter(rsqout(:,1),rsqout(:,5))
legend('R All data','R Central','R Ascending','R Descending')
xlabel('Spring Setting');
ylabel('R Value');

%% Plot of regression lines across the range of spring settings
heightplot=linspace(minheight,maxheight,(maxheight-minheight+1));
legendset={};
figure;
for q=1:15
    if any(regressout(q,:))
        pplot=regressout(q,:)';
        forceplot=polyval(pplot,heightplot);
        hold on
        plot(heightplot,forceplot)
        legendset=[legendset int2str(q)];      
    else
    end
end
 title(' Predicted Force vs. Height for Range of Spring Settings');
    ylabel('Force, N');
    xlabel('Height, mm');
    legend(legendset);
    title(legend,'Spring Setting')
    %%  Remove data for Spring Setting 15
    
    %  Inspection of the plots shows that Spring setting 15 appears
    %  inconsistent, probably because the force involved caused the MAS to
    %  tip.  As such, this data will not be used further.  
    regressout(15,:)=0;
    %%  Generate 3D plot of force, height and spring setting
    
    %  As a further visualisation of the relationship
    
    %surfplot=zeros(length(heightplot),3);
    springset=[];
    forceplotout=[];
    
    for d= 1:15
        if any(regressout(d,:))
        pplot=regressout(d,:)';
        forceplot=polyval(pplot,heightplot);
        springset= [springset;d];
        forceplotout=[forceplotout;forceplot];
        else
        end
    end
    
    % plot the 3D data
  figure;
  X=heightplot;
  Y=springset;
  Z=forceplotout;
  meshc(X,Y,Z);
  view([-80,20]);
  ylabel('Spring Setting');
  xlabel('Height, mm');
  zlabel('Force, N');
  title('Predicted Force vs. Height over range of Spring Settings')
%% Generate expression to model force as fn(height,spring setting)
    
    % Now generate a single expression to predict the force at a height
    % within the range, given the spring setting
    % First find a line to describe how forces at mid range of height  
    % vary with spring setting

forcemidrange=forceplotout(:,round(length(forceplotout)/2));
pmid=polyfit(springset,forcemidrange,1);
forcemidcalc=polyval(pmid,springset);
%Check this linear expression is a good fit for the original data
figure;
plot(springset,forcemidrange,springset,forcemidcalc)
  
%  Now use this expression for the mid-height force, and the average slope
%  from the range of spring settings, to generate a new simplified mesh of
%  spring-setting, height and force
%  First find the average slope of the force:height lines for the various
%  Spring settings:

coeffs=[];
   for d= 1:15
        if any(regressout(d,:))
            coeffs=[coeffs;regressout(d,:)];
        else
        end
   end
avslope=mean(coeffs(:,1));

%  Next calculate the force at mid-height range for all spring settings
springrange=linspace(1,15,15)';
forcemid=polyval(pmid,springrange);
%  Then extrapolate back with the average slope to find the forces at the
%  minimum height for each spring setting
forcemin=forcemid-(avslope*round(length(forceplotout)/2));
% and find the line which describes these minimum forces
pmin=polyfit(springrange,forcemin,1);
% calculate this line and add it to the plot
forcemincalc=polyval(pmin,springrange);
%  Add these minimum height forces to the plot : NOT IN REPORT
%hold on
%plot(springrange,forcemincalc)

  %legend('Actual Mid-range Force','Estimated Mid-range Force','Estimated Min Range Force');
  legend('Actual mid-range Force','Estimated mid-range Force');
  xlabel('Spring Setting');
  ylabel('Force, N');
  title('Force at mid-range height: Original model and linear estimate')

 %% Compare new model with previous results in 3D plot
 
 %  Generate a new 3D set of data using the minimum force line as the 
 %  starting point, and extrapolating forward using the average slope
  Xx=heightplot;
  Yy=springrange;
  Xxx=Xx-minheight+1;
  Zz(Yy,Xxx)=forcemincalc(Yy)+avslope*(Xx-minheight);
  
  %  Compare the original 3D data with this newly calculated data
  figure;
  C=Zz.*3;
  meshc(Xx,Yy,Zz,C)
  view([-80,20]);
  X=heightplot;
  Y=springset;
  Z=forceplotout;
  hold on
  meshc(X,Y,Z);
  ylabel('Spring Setting');
  xlabel('Height, mm');
  zlabel('Force, N');
  %legend('Original Data','Model');
  title({'Predicted Force vs. Height over range of Spring Settings','Blue plot is original data, Orange/yellow is model'})
  
%% Compare new model with previous results in 2D plot

%  Repeat with the 2D plots of force vs. height for various spring settings
%  as its easier to compare with the original

legendset={};
figure;

for q=1:15
    if any(regressout(q,:))
        pplot=regressout(q,:)';
        forceplot=polyval(pplot,heightplot);
        hold on
        plot(heightplot,forceplot)
        legendset=[legendset int2str(q)];      
    else
    end
end
     for n=1:15
        hold on
        plot(heightplot,Zz(n,:),'c--')
     end
    title(' Original vs. Predicted Force : Height Relationship for Range of Spring Settings');
    ylabel('Force, N');
    xlabel('Height, mm');
    legendset=[legendset 'Modelled']; 
    legend(legendset);
    title(legend,'Spring Setting')
%%  Plot the new model alone in 2D

% Show just the predicted forces in a 2D plot for clarity
legendset={};
figure;
for n=1:15
        hold on
        plot(heightplot,Zz(n,:))
        legendset=[legendset int2str(n)];  
     end
    title(' Predicted Force vs. Height for Range of Spring Settings');
    ylabel('Force, N');
    xlabel('Height, mm');
    legend(legendset);
    title(legend,'Spring Setting')
    
%% Compare original data, central regression line and new model

%  Repeat the earlier plots of the original data, with the new model
%  included

for cc=1:15
    if ~isempty(rangeout{1,cc})
        datacell=rangeout{1,cc};
        height=datacell(:,1);
        force=datacell(:,2);
        % fit linear regression line to data
        p=polyfit(height,force,1);
        % find predicted value
        calcforce=polyval(p,height);
        %plot results        
        figure;
        plot(height,force)
        hold on
        %plot (height,calcforce)
        
        % use the split data into ascending and descending parts
        
        heightasc=ascending{1,cc}(:,1);
        forceasc=ascending{1,cc}(:,2);
        heightdesc=descending{1,cc}(:,1);
        forcedesc=descending{1,cc}(:,2);
        % create separate regression lines for ascending and descending
        pasc=polyfit(heightasc,forceasc,1);
        calcforceasc=polyval(pasc,heightasc);
        pdesc=polyfit(heightdesc,forcedesc,1);
        calcforcedesc=polyval(pdesc,heightdesc);
        %plot those lines
        hold on
       % plot(heightasc,calcforceasc)
        hold on
        %plot(heightdesc,calcforcedesc)
        
        %  create a new set of coefficients for a central line by combining
        %  the ascending and descending coefficients
        pcentral=(pasc+pdesc)/2;
        calcforcecentral=polyval(pcentral,height);
        hold on
        plot(height,calcforcecentral)
        
        %  Now add the line generated by the new model
        hold on
        plot(heightplot,Zz(cc,:))
        
        ylabel('Force, N');
        xlabel('Height, mm');
        legend('Original Data','LS Central','New Model')
        title(['Spring setting ',int2str(cc),' Force vs. Displacement']);
        
    else 
    end
end
%% Compare original data and new model

%  And finally, plot just the original data and the new model..

for cc=1:15
    if ~isempty(rangeout{1,cc})
        datacell=rangeout{1,cc};
        height=datacell(:,1);
        force=datacell(:,2);
            
        figure;
        plot(height,force)
                
        %  Now add the line generated by the new model
        hold on
        plot(heightplot,Zz(cc,:))
        
        ylabel('Force, N');
        xlabel('Height, mm');
        legend('Original Data','New Model')
        title(['Spring setting ',int2str(cc),' Force vs. Displacement']);
        
    else 
    end
end
%%  Analyse hysteresis data

% remove final row of hysteresis data from spring setting 15 - not reliable
hysteresis(length(hysteresis),:)=[];
%%
% collate friction force, assuming its half of total hysteresis
minheightfriction=0.5*(hysteresis(:,2)-hysteresis(:,4));
minheightcentralforce=-hysteresis(:,3);
maxheightfriction=0.5*(hysteresis(:,5)-hysteresis(:,7));
maxheightcentralforce=-hysteresis(:,6);
springsetting=hysteresis(:,1);

figure;
plot(springsetting,minheightfriction)
hold on
plot(springsetting,maxheightfriction)
plot(springsetting,minheightcentralforce)
plot(springsetting,maxheightcentralforce)

ylabel('Force N');
xlabel('Spring Setting');
legend('Friction at min height','Friction at max height','Central force at min height','Central force at max height');
title('Variation of friction and central force with height and spring setting');

frictiondata=[minheightfriction;maxheightfriction];
forcedata=[minheightcentralforce;maxheightcentralforce];

figure;
scatter(forcedata,frictiondata)
xlabel('Arm Force, N');
ylabel('Friction, N');
title('Average arm force vs. frictional force in system');
pfric=polyfit(forcedata,frictiondata,1);
% find predicted value
calcfric=polyval(pfric,forcedata);
hold on
plot(forcedata,calcfric)

 % Find R value to gauge fit
        % find residual
        fricresid=frictiondata-calcfric;
        % find R value
        SSfricresid=sum(fricresid.^2);
        SSfrictotal = (length(frictiondata)-1) * var(frictiondata);
        rsq = 1 - SSfricresid/SSfrictotal;
        pfric;
        rsqstr=sprintf('%0.2f',rsq);
        pfricstr=sprintf('%0.2f',pfric(1));
        textstr=['Rsq for regression line is ',rsqstr];
        text(10,20,textstr);
        textstr2=['Slope is ',pfricstr];
        text(10,18,textstr2);
            
        %friction=frictiondata./forcedata

%%  Conclusion:  The model

% So the new model can be described by the equation:
%
% Force = -5.778*SpringSetting +0.0073*(Height-500) -3.854   

% But as the sign convention is opposite in the forward dynamic model,
%
% Force(FDmodel) = 5.778*SpringSetting-0.0073*(Height-500)+3.854

%  Also, the correlation between hysteresis force and arm force (Rsq=0.95)
%  supports the hypothesis that the hyseresis is due to frictionin the
%  system, and that this friction may be modelled as 21% of the average
%  force for that spring setting and height.

