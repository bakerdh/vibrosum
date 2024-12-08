function tactiledippers

% experiment code to run dipper experiment using 10-channel vibrotactile device
% responses are with the 3-channel footpedal
% first checks which USB audio device the headphones are connected to
% each block lasts around 6 minutes, and participant progress is indicated
% after each block
% DHB 19/10/23

close all;

subj = menu('Select participant','SW','DHB','CEJP','Test');

E.exptpath = strcat(pwd, '/');
KbName('UnifyKeyNames');   % probably unnecessary, but standardises keyboard characters across operating systems
ListenChar(-1);             % turn off keyboard echoing to stop script filling up with characters from the footpedal
rng('shuffle');         % reseed the random number generator from the clock

ST.duration = 0.5;
ST.ISI = 0.4;
ST.F1 = 26;
ST.beepfreq = 440;
ST.postresponseduration = 0.2;

E.pedcontrastsdB = -12:6:30;
E.pedcontrastsC = (10.^(E.pedcontrastsdB./20))./100;
E.pedcontrastsC(1) = 0;

nSCs = [4 8 8 8 8 8 8 8];

subnamelist = {'P1','P2','P3','XX'};

if ~exist('Results','dir')
    mkdir('Results');
end
if ~exist(strcat(E.exptpath,'Results/', subnamelist{subj}, '/'),'dir')
    mkdir(strcat(E.exptpath,'Results/', subnamelist{subj}, '/'));
end

% here generate or load in a results file to keep track of progress
if exist(strcat('Results/',subnamelist{subj},'conditionorder.mat'),'file')
    load(strcat('Results/',subnamelist{subj},'conditionorder.mat'));
    disp('Reloaded participant settings');
else
    condorder = [randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC)) randperm(length(E.pedcontrastsC))];
    currentrep = 0;
    
    save(strcat('Results/',subnamelist{subj},'conditionorder.mat'),'condorder','currentrep');
    disp('Saved participant settings');
end

% first find all available audio devices, and get ID cnumbers for the
% 6-channel USB audio devices
devices = PsychPortAudio('GetDevices');
idlist = [];
samplerates = [];
for d = 1:length(devices)
    if (devices(d).NrOutputChannels==6)
        idlist(end+1) = devices(d).DeviceIndex;
        samplerates(end+1) = devices(d).DefaultSampleRate;
    end
end
if length(idlist)<1
    disp('No multichannel audio devices detected!');
end
if length(idlist)==1
    disp('Only one multichannel audio device detected!');
end

InitializePsychSound(1);
device1 = PsychPortAudio('Open',idlist(1),1,3,44100,6);
device2 = PsychPortAudio('Open',idlist(2),1,3,44100,6);


% set up onscreen window for feedback
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
oldSuppressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
screens=Screen('Screens');
screenNumber=min(screens);
rect = [1 1 1024 1024];
[w, winRect] = Screen('OpenWindow', screenNumber, 0);  %, rect);
Screen('FillRect', w, 0);
Screen('Flip', w);
[width, height] = Screen('WindowSize', w);

HideCursor;

feedbackboxwidth = 200;
r1 = [1 1 feedbackboxwidth feedbackboxwidth];
destRect = CenterRectOnPoint(r1, width*0.5, height*0.5);


% next play the wav files telling the subject to press either the left or
% right foot pedal, depending on what they hear through the headphones

l = audioread('LeftPedal2.wav');
r = audioread('RightPedal2.wav');

l = repmat(0.25.*l',[6 1]);
r = repmat(0.25.*r',[6 1]);

PsychPortAudio('FillBuffer', device1, l);
PsychPortAudio('FillBuffer', device2, r);

exitcode = 0;

while exitcode==0
    
    onset = GetSecs + 0.2;
    PsychPortAudio('Start', device1, [], onset);
    PsychPortAudio('Start', device2, [], onset);
    
    
    while exitcode==0 && (GetSecs<onset+5)
        
        [x,y,buttons] = GetMouse;
        [keyIsDown, secs, keyCode] = KbCheck;
        
        if keyCode(KbName('a'))
            resp = 1;
            exitcode = 1;
        elseif keyCode(KbName('c'))
            resp = 2;
            exitcode = 1;
        elseif keyCode(KbName('Escape'))
            exitcode = 1;
            resp = 0;
        end
        
    end
    
end

PsychPortAudio('Close', device1);
PsychPortAudio('Close', device2);

if resp>0
    device1handle = idlist(3-resp);
    device2handle = idlist(resp);
    device1samplerate = samplerates(3-resp);
    device2samplerate = samplerates(resp);
    
    trialbeep = MakeBeep(ST.beepfreq,ST.duration,device2samplerate).*0.25;
    stimwaveform = MakeBeep(ST.F1,ST.duration,device2samplerate);
    
    pedlevel = condorder(currentrep+1);
    
    SC = getstaircasestruct(nSCs(pedlevel),100*E.pedcontrastsC(pedlevel));
    R.ntrials = zeros(SC.ncases,length(SC.levels));
    R.ncorrect = R.ntrials;
    alltrials.subj = subnamelist{subj};
    
    
    WaitSecs(5);
    
    
    trialcounter = 0;
    while sum(SC.finished)<SC.ncases		% do trials if staircases are incomplete
        
        trialcounter = trialcounter + 1;
        testinterval = ceil(rand*2);      % test in interval 1 or 2
        alltrials.testinterval(trialcounter) = testinterval;
        
        currenttrial = 0;
        while currenttrial==0                   % select current condition
            trialcond = ceil(rand*SC.ncases);
            if SC.finished(trialcond)==0
                currenttrial = trialcond;
            end
        end
        
        alltrials.staircase(trialcounter) = currenttrial;
        alltrials.pedcontrast(trialcounter) = 100*E.pedcontrastsC(pedlevel);
        alltrials.targetcontrast(trialcounter) = 100*SC.levelsM(currenttrial,SC.currentno(currenttrial));
        
        % generate stimuli for this trial
        switch currenttrial
            case 1      % monocular left
                TLstim = stimwaveform.*(E.pedcontrastsC(pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, set A
                TRstim = stimwaveform.*0;                                           % target interval, set B
                PLstim = stimwaveform.*E.pedcontrastsC(pedlevel);              % null interval, set A
                PRstim = stimwaveform.*0;                                           % null interval, set B
            case 2      % monocular right
                TLstim = stimwaveform.*0;                                   % target interval, set B
                TRstim = stimwaveform.*(E.pedcontrastsC(pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, set A
                PLstim = stimwaveform.*0;                                   % null interval, set B
                PRstim = stimwaveform.*E.pedcontrastsC(pedlevel);      % null interval, set A
            case {3,4}  % binocular
                TLstim = stimwaveform.*(E.pedcontrastsC(pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, set A
                TRstim = stimwaveform.*(E.pedcontrastsC(pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, set A
                PLstim = stimwaveform.*E.pedcontrastsC(pedlevel);      % null interval, set A
                PRstim = stimwaveform.*E.pedcontrastsC(pedlevel);      % null interval, set A
            case 5      % hbin left
                TLstim = stimwaveform.*(E.pedcontrastsC(pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));      % target interval, set A
                TRstim = stimwaveform.*E.pedcontrastsC(pedlevel);                            % target interval, set B
                PLstim = stimwaveform.*E.pedcontrastsC(pedlevel);                            % null interval, set A
                PRstim = stimwaveform.*E.pedcontrastsC(pedlevel);                            % null interval, set B
            case 6      % hbin right
                TLstim = stimwaveform.*E.pedcontrastsC(pedlevel);                            % target interval, set A
                TRstim = stimwaveform.*(E.pedcontrastsC(pedlevel)+SC.levelsM(currenttrial,SC.currentno(currenttrial)));       % target interval, set B
                PLstim = stimwaveform.*E.pedcontrastsC(pedlevel);                            % null interval, set A
                PRstim = stimwaveform.*E.pedcontrastsC(pedlevel);                            % null interval, set B
            case 7      % dich left
                TLstim = stimwaveform.*SC.levelsM(currenttrial,SC.currentno(currenttrial));      % target interval, set A
                TRstim = stimwaveform.*E.pedcontrastsC(pedlevel);                           % target interval, set B
                PLstim = stimwaveform.*0;                                                        % null interval, set A
                PRstim = stimwaveform.*E.pedcontrastsC(pedlevel);                           % null interval, set B
            case 8      % dich right
                TLstim = stimwaveform.*E.pedcontrastsC(pedlevel);                            % target interval, set B
                TRstim = stimwaveform.*SC.levelsM(currenttrial,SC.currentno(currenttrial));       % target interval, set A
                PLstim = stimwaveform.*E.pedcontrastsC(pedlevel);                            % null interval, set B
                PRstim = stimwaveform.*0;                                                         % null interval, set A
        end
        
        Tmat1 = zeros(6,length(stimwaveform));
        Tmat1(1,:) = TLstim;
        Tmat1(3,:) = TLstim;
        Tmat1(5,:) = TLstim;
        Tmat1(2,:) = TRstim;
        Tmat1(4,:) = TRstim;
        Tmat1(6,:) = TRstim;
        Tmat2 = Tmat1;
        Tmat2(5,:) = trialbeep;
        Tmat2(6,:) = trialbeep;
        
        Pmat1 = zeros(6,length(stimwaveform));
        Pmat1(1,:) = PLstim;
        Pmat1(3,:) = PLstim;
        Pmat1(5,:) = PLstim;
        Pmat1(2,:) = PRstim;
        Pmat1(4,:) = PRstim;
        Pmat1(6,:) = PRstim;
        Pmat2 = Pmat1;
        Pmat2(5,:) = trialbeep;
        Pmat2(6,:) = trialbeep;
        
        % make buffers to store stimuli in
        int1a = PsychPortAudio('Open',device1handle,1,3,device1samplerate,6);
        int1b = PsychPortAudio('Open',device2handle,1,3,device2samplerate,6);
        int2a = PsychPortAudio('Open',device1handle,1,3,device1samplerate,6);
        int2b = PsychPortAudio('Open',device2handle,1,3,device2samplerate,6);
        
        % load waveforms into buffers
        if testinterval==1
            PsychPortAudio('FillBuffer', int1a, Tmat1);
            PsychPortAudio('FillBuffer', int1b, Tmat2);
            PsychPortAudio('FillBuffer', int2a, Pmat1);
            PsychPortAudio('FillBuffer', int2b, Pmat2);
        else
            PsychPortAudio('FillBuffer', int1a, Pmat1);
            PsychPortAudio('FillBuffer', int1b, Pmat2);
            PsychPortAudio('FillBuffer', int2a, Tmat1);
            PsychPortAudio('FillBuffer', int2b, Tmat2);
        end
        
        onset = GetSecs + 0.2;
        PsychPortAudio('Start', int1a, [], onset);
        PsychPortAudio('Start', int1b, [], onset);
        
        WaitSecs(ST.duration + ST.ISI - 0.2);
        
        onset = GetSecs + 0.2;
        PsychPortAudio('Start', int2a, [], onset);
        PsychPortAudio('Start', int2b, [], onset);
        
        WaitSecs(ST.duration);
        
        prdstart = GetSecs;
        exitcode = 0;
        breakcode = 0;
        while exitcode==0
            
            [x,y,buttons] = GetMouse;
            [keyIsDown, secs, keyCode] = KbCheck;
            
            if keyCode(KbName('a'))
                resp = 1;
                exitcode = 1;
                responsetime = GetSecs - prdstart;
            elseif keyCode(KbName('c'))
                resp = 2;
                exitcode = 1;
                responsetime = GetSecs - prdstart;
            elseif sum(buttons)
                exitcode = 1;
                breakcode = 1;
                SC.finished(:) = 1;
                responsetime = 0;
                resp = 0;
            end
            
        end
        
        Screen('FillRect', w, 0);
        Screen('Flip',w);
        
        R.ntrials(currenttrial,SC.currentno(currenttrial)) = R.ntrials(currenttrial,SC.currentno(currenttrial)) + 1;		% update results structure
        if resp==testinterval
            alltrials.iscorrect(trialcounter) = 1;
            R.ncorrect(currenttrial,SC.currentno(currenttrial)) = R.ncorrect(currenttrial,SC.currentno(currenttrial)) + 1;
            SC = dostaircase(SC, currenttrial, 1);      % correct
            %             PsychPortAudio('Start', PPA.gb);
            disp('Correct');
        else
            alltrials.iscorrect(trialcounter) = 0;
            SC = dostaircase(SC, currenttrial, 0);      % incorrect
            %             PsychPortAudio('Start', PPA.bb);
            disp('Incorrect');
        end
        alltrials.responsetime(trialcounter) = responsetime*1000;
        
        PsychPortAudio('Close', int1a);
        PsychPortAudio('Close', int1b);
        PsychPortAudio('Close', int2a);
        PsychPortAudio('Close', int2b);
        
        Screen('FillRect', w, 0);
        if resp==testinterval
            Screen('FillRect', w, [0 255 0], destRect);
        else
            Screen('FillRect', w, [255 0 0], destRect);
        end
        Screen('Flip',w);
        WaitSecs(ST.postresponseduration);
        
        Screen('FillRect', w, 0);
        Screen('Flip',w);
        
    end
    
    ListenChar(1);  % turn keyboard echoing back on
    
    if ~breakcode
        currentrep = currentrep + 1;
        save(strcat('Results/',subnamelist{subj},'conditionorder.mat'),'condorder','currentrep');
        R.levels = SC.levels;
        save(strcat(E.exptpath,'Results/', subnamelist{subj}, '/P',num2str(pedlevel),'rep',num2str(ceil(currentrep/length(E.pedcontrastsdB))),'.mat'), 'R','alltrials');
    end
    
    fid = fopen(strcat(E.exptpath,'Results/', subnamelist{subj}, '/P',num2str(pedlevel),'rep',num2str(ceil(currentrep/length(E.pedcontrastsdB))),'.csv'),'w');
    fprintf(fid,'Subject,Condition,PedestalContrast,TargetContrast,TargetInterval,IsCorrect,ResponseTime\n');
    for s = 1:length(alltrials.testinterval)
        outputvect = [alltrials.staircase(s), alltrials.pedcontrast(s), alltrials.targetcontrast(s), alltrials.testinterval(s), alltrials.iscorrect(s), alltrials.responsetime(s)];
        fprintf(fid,strcat(alltrials.subj, ',%2.0f,%2.3f,%2.3f,%2.0f,%2.0f,%2.3f\n'),outputvect);
    end
    fclose(fid);
    
end

Screen('CloseAll');
ShowCursor;

end
%--------------------------------------------------------------------------------------------------
function SC = getstaircasestruct(ncases,pedcontrast)

% now with separate staircase levels for different conditions
% dichoptic can get up to 100% (40dB)
% other conditions up to the closest dB value to 100-pedestal

maxval = floor(20*log10(100-pedcontrast));  % maximum dB value
minval = maxval - 57;
SC.pedlevelC = pedcontrast;
SC.ncases = ncases;                % staircase pairs
SC.stepsize = 3;		% staircase step size (dB)
SC.downrule = 3;		% decreases contrast after x correct responses
SC.uprule = 1;			% increases contrast after x incorrect responses
SC.maxtrials = 70;
SC.maxreversals = 12;   %
SC.minlevel = [minval minval minval minval minval minval -17 -17];
SC.maxlevel = [maxval maxval maxval maxval maxval maxval 40 40];
for cond = 1:SC.ncases
    SC.levels(cond,:) = SC.minlevel(cond):SC.stepsize:SC.maxlevel(cond);	   % possible test contrast values
end
SC.levelsM = 10.^(SC.levels./20)/100;                  % contrast levels in michelson
SC.startpoint(1:ncases) = 14;   % starting point in 'staircase units'
for cond = 1:SC.ncases
    SC.startlev(cond) = SC.levels(cond,SC.startpoint(cond));  % starting point in dB
end
SC.ntrials(1:SC.ncases) = 0;			% reset staircase variables for this block
SC.nreversals(1:SC.ncases) = 0;
SC.finished(1:SC.ncases) = 0;
SC.nright(1:SC.ncases) = 0;
SC.nwrong(1:SC.ncases) = 0;
SC.currentno(1:SC.ncases) = SC.startpoint;		% starting point for the staircase (staircase units)
SC.lastdir(1:SC.ncases) = 0;

end
%--------------------------------------------------------------------------
function SC = dostaircase(SC, currenttrial, response)

% this adjusts the staircase structure (SC) so that it is correct for the next trial
% current trial gives the condition (here 1-4) of the trial which has just been run
% response is either 0 (incorrect) or 1 (correct)

SC.ntrials(currenttrial) = SC.ntrials(currenttrial) + 1;		% trial counter
jump = 1;
if SC.nreversals(currenttrial)<2    % after the first trial nreverals = 1, so setting this to 2 means we really start collecting on the 1st reversal
    jump = 2;                                                   % goes in bigger steps if we're before the first reversal
    SC.ntrials(currenttrial) = SC.ntrials(currenttrial) - 1;    % don't count trials before the first reversal
end

if response==0					% add to the numbers of right or wrong responses at this level
    SC.nwrong(currenttrial) = SC.nwrong(currenttrial) + 1;
    % SC.nright(currenttrial) = 0;
else
    SC.nright(currenttrial) = SC.nright(currenttrial) + 1;
    % SC.nwrong(currenttrial) = 0;
end

thisdir = SC.lastdir(currenttrial);
if SC.nwrong(currenttrial)==SC.uprule			% need to increment
    SC.currentno(currenttrial) = SC.currentno(currenttrial) + jump;
    SC.nwrong(currenttrial) = 0;                    % reset to 0
    thisdir = 1;									% going up
elseif SC.nright(currenttrial)==SC.downrule		% need to decrement
    SC.currentno(currenttrial) = SC.currentno(currenttrial) - jump;
    SC.nright(currenttrial) = 0;                    % reset to 0
    thisdir = -1;									% going down
end

if thisdir~=SC.lastdir(currenttrial)		% have we changed direction?
    SC.nreversals(currenttrial) = SC.nreversals(currenttrial) + 1;
    SC.lastdir(currenttrial) = thisdir;
end

if SC.currentno(currenttrial)>length(SC.levels)			% check to see if we've gone too high
    SC.currentno(currenttrial) = length(SC.levels);
elseif SC.currentno(currenttrial)<1						% or too low
    SC.currentno(currenttrial) = 1;
end

if SC.ntrials(currenttrial)>=SC.maxtrials				% has the staircase finished?
    SC.finished(currenttrial) = 1;
elseif SC.nreversals(currenttrial)>=SC.maxreversals
    SC.finished(currenttrial) = 1;
end

end
%--------------------------------------------------------------------------
