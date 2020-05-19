function [] = tdtExtract(filepath1,blockname1,savename1,polyOrder)

% Code written by Christian Pedersen
% Michael Bruchas Lab - UW
% Code adapted by Jennifer Dezet Deem and Sarah Larsen

% Extracts photometry data from TDT Data Tank and outputs time series as
% .mat file

% REQUIRES tdt2mat.m file to be in same folder (filepath)

%% 

clc   % clear command window


%%%%%%%%%%%%%%%%%%%%%%  EDIT THESE FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

tankdir = filepath1;
tankname = 'AgRP-test-190924-105627'; % name of your tank

%blockname = 'Pacap15-4-191018-120914'; % name of your file
blockname = blockname1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
%LMag is the demodulated data, may also have other timestamps etc

storenames2 = {'405A'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename);
end

for k = 1:numel(storenames2)
  storename2 = storenames2{k};
  S2{k} = tdt2mat(tankdir, tankname, blockname, storename2);
end

% Massage data and get time stamps

LMag = S{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani1 = LMag.channels==1;
chani2 = LMag.channels==2;

LMag2 = S2{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani21 = LMag2.channels==1;
chani22 = LMag2.channels==2;
% chani21 = LMag2.channels==1;
% chani22 = LMag2.channels==2;

% Get LMag data as a vector (repeat for each channel)
rawdat1 = LMag.data(chani1,:);
rawdat1 = reshape(rawdat1', [],1); % unwrap data from m x 256 array
% dat2 = LMag.data(chani21,:);
% dat2 = reshape(dat2', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts = LMag.timestamps(chani1);
t_rec_start = ts(1);

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
ts = reshape(ts',[],1);

%%%%%%%%%%%%%%%%%%

dat2 = LMag2.data(chani21,:);
dat2 = reshape(dat2', [],1); % unwrap data from m x 256 array
% dat2 = LMag.data(chani21,:);
% dat2 = reshape(dat2', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts2 = LMag2.timestamps(chani21);
t_rec_start2 = ts2(1);

ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
ts2 = reshape(ts2',[],1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subtracted 405 signal from 470 signal % raw final signal
subdat = rawdat1-dat2;
Dts=ts;

% plot raw 405 ch, raw 470 ch and subtracted signal
figure()
hold on
plot(ts,rawdat1,'b');
plot(ts2,dat2,'r');
title('both signals')
plot(ts,subdat,'g');
xlabel('time(s)')
ylabel('amplitude')
hold off


%% fit 2nd order exponential or polynomial to estimate and subtract baseline drift

%f2 = fit(Dts,subdat,'exp2');
%fitcurve= f2(Dts);

p1 = polyfit(Dts,subdat,polyOrder);
fitcurve = polyval(p1,Dts);

dataFilt = subdat - fitcurve;


%% calculate dF/F relative to first 5 minutes of baseline

Fs = round(1017.25);

normDat1 = (dataFilt - median(dataFilt(1:(300*Fs))))./abs(median(rawdat1(1:(300*Fs))));
normDat = normDat1.*100; % make a percentage

data1 = normDat;

% plot subtracted signal w/o baseline correction
figure();
plot(Dts,subdat,Dts,fitcurve);
title('Subtracted signal')
xlabel('time(s)')
ylabel('raw F')

% plot final, baseline-corrected signal
figure();
plot(Dts,data1);
title('Signal with corrected baseline')
xlabel('time(s)')
ylabel('dF/F (%)')
ylim([-20 20])

%% Save file as .mat file with specified filename

save(savename1,'data1','Dts');



end





