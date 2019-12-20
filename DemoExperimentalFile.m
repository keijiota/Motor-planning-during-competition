clear; clc;
% Screen('Preference', 'SkipSyncTests', 1);

%% initialization
scNumber = max(Screen('Screens'));
bgColor = [0 0 0];
[wPtr, winRect] = Screen('OpenWindow', scNumber,  bgColor , [0 0 500 500]) ;

Priority(1);
KbName('UnifyKeyNames'); % get key code
spaceKeyID = KbName('space');
Screen('BlendFunction', wPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % for anti-aliasing
Screen('SelectStereoDrawBuffer', wPtr, 0);
Screen('TextSize', wPtr, 22); 
font = {'Arial'}; 
Screen('TextFont', wPtr, char(font));

pixelpermm = mm2pixeltarget(0, winRect) - mm2pixeltarget(1, winRect);

%% define states
closeexp = 0;
trialStart = 1;
nextstate = 2;

%% main code
ellapsedTime = 0 ;
state = trialStart ;

iTrial = 0;
nTrial = 2;

while state ~= closeexp
    
    % ----time stamp----
    timeStamp = GetSecs ;
    [x y button] = GetMouse(0);
    x
    button
    
    
    %----switch state when given condition is satisfied-----------
    switch state
        case trialStart
            iTrial = iTrial+1 ;
            if iTrial <= nTrial
                state = nextstate ;
                time_nextstate = timeStamp; % get current time
            else
                state = closeexp ;
            end
            
        case nextstate
            if timeStamp - time_nextstate > 50 
                state = trialStart;
            end
    end
    
    %----- we draw stimuli below--------
    
    % Draw Trials
    Screen('DrawText',wPtr, sprintf('Trial %d/%d', iTrial, nTrial), [winRect(3)/2],[winRect(4)/2],[255 255 255]);
    
    %-----Draw all stimuli here defined above---------
    [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', wPtr);
    
    %------Press space key whenever you want to stop the program
    [keyIsDown, secs, keycode] = KbCheck;
    if keycode(spaceKeyID) == 1
        state = closeexp ;
    end
end 


%% finalization
Screen('CloseAll');
% filename = [Filename,'_All'] ;
% save(filename, 'DataAllTrial');


 