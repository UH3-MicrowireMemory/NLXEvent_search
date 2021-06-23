% NLX event search script
%
% Identify CSC file that matches with event files for a participant.
%
% Instructions: 1. Find event file in saveSessAll that matches timestamp with .txt file of
% behavioral data recording. 2. Enter that saveSess # into saveSessAll{_}.
% 
% Need to change: 1. Directory location of 'allSessionData.mat' 
%                 2. Directory location of CSC files
%                 3. Make sure Nlx2Mat functions are on path 
% 
% Output 'SessionDATE' = Timestamp for event file to compare with .txt file
% Output 'foundCSC' = CSC file name that matches the event file
%
% John A. Thompson | May 28 2021

% Add Nlx2Mat functions to path
addpath('C:\Users\darwinm\Documents\Github\NLX-Event-Viewer\NLX_IE_Code');

% Load location of 'allSessionData.mat' for patient
cd('C:\Users\darwinm\Documents\Thompson Lab\Microwire\PatientData\MW1\SessionSaves');

% Locate saveSess file that matches timestamp by plugging file #'s in {} 
timeSTAMP=saveSessAll{10}.StartTime;
timeREcord = datetime(timeSTAMP/1000000,'ConvertFrom','posixtime','TimeZone','America/Denver');
SessionDATE=timeREcord+hours(6)

%cd to CSC files for patient
cd('C:\Users\darwinm\Documents\Thompson Lab\Microwire\PatientData\MW1\Events');

% Locate unique CSC extension files (i.e, CSC1__001, __002)
cscRepoT = dir('*.ncs');
cscRepo = {cscRepoT.name};

fileEls = cellfun(@(x) strsplit(x, {'_','.'}), cscRepo, 'UniformOutput',false);

cscAllnum = cellfun(@(x) x{1}, fileEls , 'UniformOutput',false);

uniCSC = unique(cscAllnum);
useCSC = uniCSC(1);

allCSCns = transpose(cscRepo(ismember(cscAllnum,useCSC)));

% Loop through all extension
micAll = zeros(size(allCSCns));
ephInd = zeros(size(allCSCns));

for ci = 1:length(allCSCns)

   [TimestampsCSC, ~, ~, ~,... 
       ~, ~] = Nlx2MatCSC(allCSCns{ci}, [1 1 1 1 1], 1, 1, [] );

    % Create interpolation
    TimestampsINT = nlxEVinterp(TimestampsCSC);

    % Create a table with the time difference between nearest value
    % and timestamp
    offSETtime = abs(TimestampsINT - timeSTAMP);
    [micOffset , ephysInd] = min(min(offSETtime));

    micAll(ci) = micOffset;
    ephInd(ci) = ephysInd;

end

% Output CSC duplicate with the lowest value
chckTable = table(allCSCns,micAll,ephInd);
chckTable.msAll = chckTable.micAll/(1e+6);

[~ , minREC] = min(chckTable.micAll);

foundCSC = chckTable.allCSCns{minREC}

%timeDIF = chckTable.msAll(minREC);
     
function timeMAT = nlxEVinterp(timeVEC)

      timeMAT = zeros(512,size(timeVEC,2));
      timeMAT(1,:) = timeVEC;
            for ti = 1:size(timeVEC,2)

                if ti == size(timeVEC,2)
                    tmpl = linspace(timeVEC(ti),timeVEC(ti)+16000,511);
                else
                    tmpl = linspace(timeVEC(ti),timeVEC(ti + 1),511);
                end
                tmplt = transpose(tmpl);
                timeMAT(2:512,ti) = tmplt;
            end
end