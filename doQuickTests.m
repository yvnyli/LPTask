





% for my mac
Screen('Preference', 'SkipSyncTests', 1);   


 
quickTest()
function quickTest()
% Clear the workspace
close all;
clearvars;
sca;
cleanob = onCleanup(@() mycleanup);
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);
 
% Define black and white
white = WhiteIndex(screenNumber); 
black = BlackIndex(screenNumber);
grey = white / 2;
inc = white - grey;

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);
width = windowRect(3);
height = windowRect(4);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');


% make box of cards
boxOfCards = [0 0 width/5 width/10];
boxOfCards_pos = CenterRectOnPointd(boxOfCards, width/3, yCenter);
boxColor = [0 0 0];
Screen('FillRect', window, boxColor, boxOfCards_pos);
% label it as box of cards
[boxCenterx,boxCentery] = RectCenter(boxOfCards_pos);
DrawFormattedText(window, 'BOX OF CARDS', 'center', 'center', ...
  [0.9,0.9,0.9], [], [], [], [], [], boxOfCards_pos);


% What's in the box
dummyBox = [0 0 10 10];
textCenter_pos = CenterRectOnPointd(dummyBox,boxCenterx,boxCentery - height/4);
DrawFormattedText(window,['There are 100 picture cards in this box. ',...
  'Some cards have a light blue iris on them. Other cards have a dark purple iris.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],textCenter_pos);


% Press SPACEBAR to continue
instruCenter_pos = CenterRectOnPointd(dummyBox,xCenter,height*9/10);
DrawFormattedText(window,['Press ''SPACEBAR'' to continue.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],instruCenter_pos);

Screen('Flip', window);
pause(0.2)
[~,keycode,~] = KbWait;
WaitSecs(0.1)
while ~keycode(KbName('space'))
  [~,keycode,~] = KbWait;
  WaitSecs(0.1)
end

% show example cards

iris1Im = imread('./iris1.jpg');
iris2Im = imread('./iris2.jpg');
% Make the image into a texture
iris1 = Screen('MakeTexture', window, iris1Im);
iris2 = Screen('MakeTexture', window, iris2Im);
[ih,iw] = size(iris1);
imbox = [0 0 height/5/ih*iw height/5];
imbox_pos = CenterRectOnPointd(imbox,width*2/3,height/3);


% light blue
Screen('DrawTextures', window, iris1, [], imbox_pos);
DrawFormattedText(window,['A card with a light blue iris looks like this.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],textCenter_pos);
% make box of cards
Screen('FillRect', window, boxColor, boxOfCards_pos);
DrawFormattedText(window, 'BOX OF CARDS', 'center', 'center', ...
  [0.9,0.9,0.9], [], [], [], [], [], boxOfCards_pos);
% Press SPACEBAR to continue
DrawFormattedText(window,['Press ''SPACEBAR'' to continue.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],instruCenter_pos);

Screen('Flip', window);
pause(0.2)
[~,keycode,~] = KbWait;
WaitSecs(0.1)
while ~keycode(KbName('space'))
  [~,keycode,~] = KbWait;
  WaitSecs(0.1)
end

% dark purple
Screen('DrawTextures', window, iris2, [], imbox_pos);
DrawFormattedText(window,['A card with a dark purple iris looks like this.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],textCenter_pos);
% make box of cards
Screen('FillRect', window, boxColor, boxOfCards_pos);
DrawFormattedText(window, 'BOX OF CARDS', 'center', 'center', ...
  [0.9,0.9,0.9], [], [], [], [], [], boxOfCards_pos);
% Press SPACEBAR to continue
DrawFormattedText(window,['Press ''SPACEBAR'' to continue.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],instruCenter_pos);

Screen('Flip', window);
pause(0.2)
[~,keycode,~] = KbWait;
WaitSecs(0.1)
while ~keycode(KbName('space'))
  [~,keycode,~] = KbWait;
  WaitSecs(0.1)
end


% slider
pslider = [0 0 width/3 3];
pslider_pos = CenterRectOnPointd(pslider, xCenter, height*3/4);
sliderColor = [1 0 0];
tickRect = [0 0 5 20];
tickvalues = linspace(pslider_pos(1),pslider_pos(3),101);
tickind = 50;

centeredTick = CenterRectOnPointd(tickRect,tickvalues(tickind+1),pslider_pos(2)+3*0.5);
Screen('FillRect', window, sliderColor, pslider_pos);
Screen('FillRect',window,[1,1,1],centeredTick);
responseCenter_pos = CenterRectOnPointd(dummyBox,xCenter,pslider_pos(4)-height/10);
DrawFormattedText(window, ['between ?? and ??'], 'center',...
  'center', [0.1,0.1,0.1],[],[],[],[],[],responseCenter_pos);
DrawFormattedText(window,['Out of the 100 cards, how many do you think have a',...
  ' light blue iris on them?\n Indicate your guess as a range with the slider.' ],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],textCenter_pos);
% make box of cards
Screen('FillRect', window, boxColor, boxOfCards_pos);
DrawFormattedText(window, 'BOX OF CARDS', 'center', 'center', ...
  [0.9,0.9,0.9], [], [], [], [], [], boxOfCards_pos);
% Press SPACEBAR to continue
DrawFormattedText(window,['Move the slider with LEFTARROW and RIGHTARROW.\n',...
  'Press ''N'' to confirm the lower bound.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],instruCenter_pos);  Screen('Flip', window);
WaitSecs(0.08);
% prior probability
while true  
  WaitSecs(0.08);
  [~,keycode,~] = KbWait;
  WaitSecs(0.1)
  if keycode(KbName('escape')) 
    return
  end
  if keycode(KbName('N'))
    WaitSecs(0.1)
    break
  end
  if keycode(KbName('rightarrow'))
    tickind = tickind + 1;
  elseif keycode(KbName('leftarrow'))
    tickind = tickind - 1;
  end
  tickind = min(tickind,100); tickind = max(tickind,0);
  centeredTick = CenterRectOnPointd(tickRect,tickvalues(tickind+1),pslider_pos(2)+3*0.5);
  Screen('FillRect', window, sliderColor, pslider_pos);
  Screen('FillRect',window,[1,1,1],centeredTick);
  DrawFormattedText(window, ['between ',num2str(tickind),' and ??'], 'center',...
    'center', [0.1,0.1,0.1],[],[],[],[],[],responseCenter_pos);
  DrawFormattedText(window,['Out of the 100 cards in the box, ',...
    'how many do you think have a',...
    ' light blue iris on them?\nIndicate your guess as a range with the slider.' ],...
    'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],textCenter_pos);
  % make box of cards
  Screen('FillRect', window, boxColor, boxOfCards_pos);
  DrawFormattedText(window, 'BOX OF CARDS', 'center', 'center', ...
    [0.9,0.9,0.9], [], [], [], [], [], boxOfCards_pos);
  % Press SPACEBAR to continue
  DrawFormattedText(window,['Move the slider with LEFTARROW and RIGHTARROW.\n',...
    'Press ''N'' to confirm the lower bound.'],...
    'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],instruCenter_pos);  Screen('Flip', window);
end

lb = tickind;


centeredTick = CenterRectOnPointd(tickRect,tickvalues(tickind+1),pslider_pos(2)+3*0.5);
Screen('FillRect', window, sliderColor, pslider_pos);
Screen('FillRect',window,[1,1,1],centeredTick);
DrawFormattedText(window, ['between ',num2str(lb),' and ??'], 'center',...
  'center', [0.1,0.1,0.1],[],[],[],[],[],responseCenter_pos);
DrawFormattedText(window,['Out of the 100 cards in the box, ',...
    'how many do you think have a',...
    ' light blue iris on them?\nIndicate your guess as a range with the slider.' ],...
    'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],textCenter_pos);
% make box of cards
Screen('FillRect', window, boxColor, boxOfCards_pos);
DrawFormattedText(window, 'BOX OF CARDS', 'center', 'center', ...
  [0.9,0.9,0.9], [], [], [], [], [], boxOfCards_pos);
% Press SPACEBAR to continue
DrawFormattedText(window,['Move the slider with LEFTARROW and RIGHTARROW.\n',...
  'Press ''SPACEBAR'' to confirm the upper bound.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],instruCenter_pos);  Screen('Flip', window);
WaitSecs(0.08);
while true
  WaitSecs(0.08);
  [~,keycode,~] = KbWait;
  WaitSecs(0.1);
  if keycode(KbName('escape')) 
    return
  end
  if keycode(KbName('space'))
    WaitSecs(0.1)
    break
  end
  if keycode(KbName('rightarrow'))
    tickind = tickind + 1;
  elseif keycode(KbName('leftarrow'))
    tickind = tickind - 1;
  end
  tickind = min(tickind,100); tickind = max(tickind,0);
  centeredTick = CenterRectOnPointd(tickRect,tickvalues(tickind+1),pslider_pos(2)+3*0.5);
  Screen('FillRect', window, sliderColor, pslider_pos);
  Screen('FillRect',window,[1,1,1],centeredTick);
  DrawFormattedText(window, ['between ',num2str(lb),' and ',num2str(tickind)], 'center',...
    'center', [0.1,0.1,0.1],[],[],[],[],[],responseCenter_pos);
  DrawFormattedText(window,['Out of the 100 cards, how many do you think have a',...
    ' light blue iris on them?\n Indicate your guess as a range with the slider.' ],...
    'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],textCenter_pos);
  % make box of cards
  Screen('FillRect', window, boxColor, boxOfCards_pos);
  DrawFormattedText(window, 'BOX OF CARDS', 'center', 'center', ...
    [0.9,0.9,0.9], [], [], [], [], [], boxOfCards_pos);
  % Press SPACEBAR to continue
  DrawFormattedText(window,['Move the slider with LEFTARROW and RIGHTARROW.\n',...
    'Press ''SPACEBAR'' to confirm the upper bound.'],...
    'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],instruCenter_pos);  
  Screen('Flip', window);
end

ub = tickind;




% Flip to the screen. This command basically draws all of our previous
% commands onto the screen. See later demos in the animation section on more
% timing details. And how to demos in this section on how to draw multiple
% rects at once.
% For help see: Screen Flip?
% Screen('Flip', window);

% Now we have drawn to the screen we wait for a keyboard button press (any
% key) to terminate the demo.
% For help see: help KbStrokeWait
% KbStrokeWait;

end
function mycleanup()
sca;
end