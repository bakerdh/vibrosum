function tactileEEG

% experiment code to run SSSEP experiment using 10-channel vibrotactile device
% responses are with the 3-channel footpedal
% first checks which USB audio device the headphones are connected to
% each block lasts around 7 minutes, and participant progress is indicated
% after each block
% DHB 8/2/24

close all;

subinfo = inputdlg({'Please enter participant ID:'});
subj = subinfo{1};

E.exptpath = strcat(pwd, '/');
KbName('UnifyKeyNames');   % probably unnecessary, but standardises keyboard characters across operating systems
ListenChar(-1);             % turn off keyboard echoing to stop script filling up with characters from the footpedal
rng('shuffle');         % reseed the random number generator from the clock

ST.duration = 11;
ST.ISI = 3;
ST.F1 = 26;
ST.F2 = 23;
ST.beepfreq = 440;

E.stimcontrastsdB = 12:6:36;
E.stimcontrastsC = (10.^(E.stimcontrastsdB./20))./100;
E.nconds = 12;
E.ntrialsperblock = 30;
E.nblocks = 4;
%E.ntrials = 120;
%restPeriod = 30; % Number of trials before a rest
dotriggers = 1;

if ~exist('Results','dir')
    mkdir('Results');
end

% here generate or load in a results file to keep track of progress
if exist(strcat('Results/',subj,'conditionorder.mat'),'file')
    load(strcat('Results/',subj,'conditionorder.mat'));
    disp('Reloaded participant settings');
else
    
    condlist = sort(repmat(1:12,1,5)); % replicate the condition 1:12 five times
    contrastlist = repmat(1:5,1,12); %replicate the constrast level 1:5 12 times
    
    condorder = [randperm(length(condlist)) randperm(length(condlist))]; %generates a random order of numbers from 1 to 60, repeat twice. total trials=120
    currenttrial = 0;
    
    save(strcat('Results/',subj,'conditionorder.mat'),'condorder','currenttrial','condlist','contrastlist');
    disp('Saved participant settings');
end

if dotriggers   % set up trigger device
    triggerbox = serial('/dev/tty.usbserial-BBTKUSBTTL');
    set(triggerbox,'BaudRate',115200,'DataBits',8,'Parity','none','StopBits',1,'FlowControl','none','Terminator','')
    fopen(triggerbox)
    fprintf(triggerbox,'RR')
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
    
    onset = GetSecs + 0.2;  % avoid delay, make sure the two stimui be presented at 'onset' time
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
    stimwaveform1 = MakeBeep(ST.F1,ST.duration,device2samplerate);
    stimwaveform2 = MakeBeep(ST.F2,ST.duration,device2samplerate);
    
    
    % pre-generate all stimuli and load into sound card memory
    
    for n = 1:E.nconds
        for c = 1:length(E.stimcontrastsdB)
            
            switch n
                case 1      % monocular left
                    PLstim = stimwaveform1.*E.stimcontrastsC(c);              % F1 signal, set A
                    PRstim = stimwaveform1.*0;                                % null interval, set B
                case 2      % monocular right
                    PLstim = stimwaveform1.*0;                                % null interval, set B
                    PRstim = stimwaveform1.*E.stimcontrastsC(c);              % F1 signal, set A
                case {3,4}  % binocular
                    PLstim = stimwaveform1.*E.stimcontrastsC(c);              % F1 signal, set A
                    PRstim = stimwaveform1.*E.stimcontrastsC(c);              % F1 signal, set A
                case 5      % dichoptic left
                    PLstim = stimwaveform1.*E.stimcontrastsC(c);              % F1 signal, set A
                    PRstim = stimwaveform1.*E.stimcontrastsC(4);              % F1 maskl, set B
                case 6      % dichoptic right
                    PLstim = stimwaveform1.*E.stimcontrastsC(4);                            % F1 maskl, set A
                    PRstim = stimwaveform1.*E.stimcontrastsC(c);                            % F1 signal, set B
                case 7      % Cross-monocular left
                    PLstim = stimwaveform2.*E.stimcontrastsC(c);              % F2 signal, set A
                    PRstim = stimwaveform2.*0;                                % null interval, set B
                case 8      % Cross-monocular right
                    PLstim = stimwaveform2.*0;                         % null interval, set B
                    PRstim = stimwaveform2.*E.stimcontrastsC(c);       % F2 signal, , set A
                case 9  %Cross- binocular
                    PLstim = stimwaveform1.*E.stimcontrastsC(c);      % F1 signal, set A
                    PRstim = stimwaveform2.*E.stimcontrastsC(c);      % F2 signal, set A
                case 10  % Cross-binocular
                    PLstim = stimwaveform2.*E.stimcontrastsC(c);      % F2 signal, set A
                    PRstim = stimwaveform1.*E.stimcontrastsC(c);      % F1 signal, set A
                case 11      % Cross--dichoptic-left
                    PLstim = stimwaveform1.*E.stimcontrastsC(c);                            % F1 signal, set A
                    PRstim = stimwaveform2.*E.stimcontrastsC(4);                            % F2 mask, set B
                case 12      % Cross-dichoptic right
                    PLstim = stimwaveform2.*E.stimcontrastsC(4);                            %  F2 mask, set A
                    PRstim = stimwaveform1.*E.stimcontrastsC(c);                            %  F1 signal, set B                                                  
            end
            
            % generate stimuli for this trial
            
            Pmat1 = zeros(6,length(stimwaveform1)); %generate vinration to each finger
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
            int1a(n,c) = PsychPortAudio('Open',device1handle,1,3,device1samplerate,6); %index the handle by the 12 conditions and 5 contrast levels
            int1b(n,c) = PsychPortAudio('Open',device2handle,1,3,device2samplerate,6);
            
            PsychPortAudio('FillBuffer', int1a(n,c), Pmat1);
            PsychPortAudio('FillBuffer', int1b(n,c), Pmat2);
            
        end
    end
    
    alltrials.subj = subj;
    
    for block = 1:E.nblocks
     
        disp(strcat('Ready to begin block ',num2str(block)));
        disp('Start recording, then press left or right pedal to begin');
        
        exitcode = 0;
        while ~exitcode
        [x,y,buttons] = GetMouse;
        [keyIsDown, secs, keyCode] = KbCheck;
        
        if keyIsDown
            exitcode = 1;
        end
        if sum(buttons)
            exitcode = 1;
        end
        WaitSecs(0.2);
        end
        
        
    WaitSecs(5); % before participants start,  5s to preparation
    
    
    trialcounter = 0;
    while trialcounter < E.ntrialsperblock
    %for i = 1: length(condorder) 
        trialcounter = trialcounter + 1;
        
        currenttrial = currenttrial + 1;
        condition = condlist(condorder(currenttrial)); 
        level = contrastlist(condorder(currenttrial));  
        
        alltrials.condition(trialcounter) = condition;
        alltrials.stimcontrast(trialcounter) = 100*E.stimcontrastsdB(level); 
        
        onset = GetSecs + 0.2; % avoid delay, make sure the two stimui be presented at 'onset' time
        PsychPortAudio('Start', int1a(condition,level), [], onset);
        PsychPortAudio('Start', int1b(condition,level), [], onset);
        
        WaitSecs(0.2);
        %disp(['Trial ', num2str(trialcounter), ' of ', num2str(E.ntrials)]);
        disp(['Trial ', num2str(trialcounter), ' of ', num2str(E.ntrialsperblock)]);  % Display the trial counter after every trial
                                                                                        %num2str() is a function used to convert numeric values to strings
        % % Check if it's time for a rest
        % if mod(trialcounter, restPeriod) == 0
        %      disp('Take a break! Press ''p'' to continue');
        %      while 1
        %         [~,~,kc]=KbCheck;
        %         if kc(KbName('p'))
        %             break;
        %         end
        %     end
        % end


        if dotriggers                              %send trigger in different conditions and constrast levels
            fprintf(triggerbox,'%02X',condition*10 + level);        %triggers look like: in condition1 it's: 11; 12; 13; 14; 15
        end
        
        
        alltrials.onsettime(trialcounter) = onset;
        
        WaitSecs(0.5);
        
        if dotriggers
            fprintf(triggerbox,'%02X',0);  %set the trigger back to zero
        end
        
        WaitSecs(ST.duration-0.7);
        
        [x,y,buttons] = GetMouse;
        [keyIsDown, secs, keyCode] = KbCheck;
        
        if keyIsDown
            trialcounter  = 10000;
        end
        if sum(buttons)
            trialcounter  = 10000;
        end
        
        WaitSecs(ST.ISI-0.2);
        
    end
    
    
    if trialcounter < 10000
        save(strcat('Results/',subj,'conditionorder.mat'),'condorder','currenttrial','condlist','contrastlist');
        
        fname = strcat(E.exptpath,'Results/', subj, 'trialrecord.csv');
        if exist(fname,'file')
            fid = fopen(fname,'a');
        else
            fid = fopen(fname,'w'); 
            fprintf(fid,'Subject,Condition,Contrast,OnsetTime\n'); 
        end
        
        for s = 1:length(alltrials.condition)
            outputvect = [alltrials.condition(s), alltrials.stimcontrast(s), alltrials.onsettime(s)];
            fprintf(fid,strcat(alltrials.subj, ',%2.0f,%2.3f,%2.3f\n'),outputvect);
        end
        fclose(fid);
    end
    
    end
    
end

    ListenChar(1);  % turn keyboard echoing back on

for cond = 1:12
    for level = 1:5
PsychPortAudio('Close', int1a(cond,level));
PsychPortAudio('Close', int1b(cond,level));
    end
% ShowCursor;

end
