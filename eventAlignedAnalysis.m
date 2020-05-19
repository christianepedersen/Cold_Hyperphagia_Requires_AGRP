% Plotting photometry data relative to behavioral events
% Code written by Christian Pedersen
% Code adapted by Jennifer Dezet Deem and Sarah Larsen

%%

clear all
close all
clc

%%
 
k = 1;

for sess = 1  % pick sessions to analyze and plot together


if sess == 1 %6-5 may15
photoname = 'raw_photom_1.mat';
delay1 = 210;
eventTimes = [660 1320 1980 2640 3300]+delay1;
tempTypes = [31 15 30 14 22];
end


%% pick degree of temperature to group

% choose temperature state to analyze EXAMPLE: 10, 30 or 22 for this
% experiment set
temperature = 14; % degrees celsius


if exist('tempTypes')==1

    eventTimes(tempTypes~=temperature) = [];

end    
    
%% 
  
pokeTimes1 = eventTimes;

%%

load(photoname)
experDuration = floor(max(Dts));

%%

% get photometry trace
Fs = 1017.25;
%experDuration = 1400; %seconds

C1 = zeros(experDuration,1);
%C1(dropTimes) = 1;
C1(pokeTimes) = 1;  % indicate Score times

startdelay = 0;
load(photoname)
photom1 = data1(round(Fs*startdelay)+1:round(Fs*(experDuration+startdelay)));
time = linspace(1/Fs,experDuration,experDuration*Fs);


behavior = C1;
timeBehav = linspace(0,experDuration,length(behavior));

%% Z score photometry data

photom1 = (photom1-mean(photom1))./std(photom1);

%%

figure(sess+100)

plot(timeBehav,behavior+10,'g',...
    resample(time,1,1000),resample(photom1,1,1000));

ylabel('Z scored                       Events per second')
xlabel('Time (sec)')
yticks([-10 -5 0 5 10 15 20 25 30 35 40])
yticklabels({'-10','-5','0','5','0','5','10','15','20','25','30'})

windtop = 20;
axis([0 experDuration -5 windtop])


%%

events1 = pokeTimes1;
  
timewindow = 300; % +/- seconds from event
  
%% Lick triggered averaging (of photom trace)

Fs = 1017.25; % photometry sample rate

samplewindow = round(timewindow*Fs); % samples per time period

events2 = events1(events1<(experDuration-timewindow));
events = events2(events2>timewindow);

LickTrig = zeros(length(events),2*samplewindow);

% Change lick times to random times for a control
%events = max(events).*rand(length(events),1);

for p = 1:length(events)
    
    trigindex = zeros(1,length(photom1));
    
    [~,startidx] = min(abs(time-(events(p)-timewindow)));
    trigindex(startidx:(startidx+size(LickTrig,2)-1)) = 1;
    
    LickTrig(p,:) = photom1(trigindex==1);
    
end

photoPerLick = mean(LickTrig,1);
%

sem = std(LickTrig,0,1)./sqrt(size(LickTrig,1)); % sem = std/sqrt(n)

% event triggered average
figure(sess+200)

hold on
%eb1 = errorbar(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,sem,'Color',[0.2,0.2,0.5]);
x = decimate(linspace(-timewindow,timewindow,length(photoPerLick)),500);
y = resample(photoPerLick,1,500);
eb = resample(sem,1,500);
lineProps.col{1} = 'red';
mseb(x,y,eb,lineProps,1);
%L = line([0 0],[-1 2]);
%axis([-20 20 -0.5 1])
%set(L,'Color','black')
xlabel('Peri-Event Time (sec)')
ylabel('Z score (smoothed)')
hold off

%%

itchStack{k} = LickTrig;

k = k + 1;

end

 %close all

%% combined plot

superStack = [];

for  n = 1:length(itchStack)
    superStack = cat(1,superStack,itchStack{n});
end

%LickTrig1 = superStack;

for q = 1:size(superStack,1)
    
   LickTrig1(q,:) = resample(superStack(q,:),1,300);
    
end


photoPerLick = mean(LickTrig1,1);


sem = std(LickTrig1,0,1)./sqrt(size(LickTrig1,1)); % sem = std/sqrt(n)

%% heat map (x dim: time, ydim: event, zdim: deltaF/F)
figure(6262)

%[~,idx1] = sort(mean(LickTrig1(:,60:75),2),'ascend');

%LickTrig1 = LickTrig1(idx1,:);

hold on
imagesc(linspace(-timewindow,timewindow,length(photoPerLick)),1:size(LickTrig1,1),LickTrig1)
L = line([0 0],[0 size(LickTrig1,1)+1]);
set(L,'Color','black')
xlabel('Time aligned to temp change')
ylabel('Animal')
cb = colorbar;
title(cb,'dF/F')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
xlim([-timewindow timewindow])
ylim([0 size(LickTrig1,1)+1])

hold off

%print -painters -depsc heatmap.eps

%% event triggered average
figure(6161) 

hold on
x = linspace(-timewindow,timewindow,length(photoPerLick));
y = photoPerLick;
eb = sem;
lineProps.col{1} = [0 0.5 0];
mseb(x,y,eb,lineProps,1);
%L = line([0 0],[-3 5]);
axis([-timewindow timewindow -2 2])
%set(L,'Color','black')
xlabel('Time aligned to temp change')
ylabel('dF/F')
hold off

%print -painters -depsc traceVTA.eps









