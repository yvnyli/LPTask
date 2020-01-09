% for my mac
%Screen('Preference', 'SkipSyncTests', 1);   

 
 
testCardAppearDisappear();






function LPTask()
% Clear the workspace
close all;
clearvars;
sca;
cleanob = onCleanup(@() mycleanup);
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

s = struct;

s = makeGreyBackground(s);

s = initParams(s);


s = makeBoxOfCards(s,[s.winw/6, s.winh/3]);

% s = drawBoxOfCards(s);
% s.hold_on_list = [s.hold_on_list, {'s = drawBoxOfCards(s)'}];

% What's in the box
dummyBox = [0 0 10 10];
textCenter_pos = CenterRectOnPointd(dummyBox,boxCenterx,boxCentery - height/4);
DrawFormattedText(window,['There are 100 picture cards in this box. ',...
  'Some cards have a light blue iris on them. Other cards have a dark purple iris.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],textCenter_pos);


% Press SPACEBAR to continue
instruCenter_pos = CenterRectOnPointd(dummyBox,xCenter,height*9/10);
DrawFormattedText(window,['Press ''SPACEBAR'' to continue.'],...
  'centerblock','center',[0.1,0.1,0.1],70,[],[],[],[],instruCenter_pos); %#ok<*NBRAK>

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

function s = makeGreyBackground(s)
% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);
 
% Define black and white
s.white = WhiteIndex(screenNumber); 
s.black = BlackIndex(screenNumber);
s.grey = (s.white-s.black) / 2;

% Open an on screen window
[s.window, windowRect] = PsychImaging('OpenWindow', screenNumber, s.grey);
% Maximum priority level
topPriorityLevel = MaxPriority(s.window);
Priority(topPriorityLevel);

% Query the frame duration
s.ifi = Screen('GetFlipInterval', s.window);

% Get the centre coordinate of the window
[s.xc, s.yc] = RectCenter(windowRect);
s.winw = windowRect(3);
s.winh = windowRect(4); 

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', s.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
end

function s = initParams(s)
s.card_width = s.winw/5;
s.animation_time = 0.5;
s.text_color = 0.1*s.white;
s.slider_len = s.winw/4;
s.slider_thickness = 10 / 1440 * s.winh;
end

function s = makeBoxOfCards(s, front_center)
% make box of cards
box_front = [0 0 s.winw/8 s.winw/8*0.7];
box_front_pos = CenterRectOnPointd(box_front,front_center(1),front_center(2));
box_front_color = 0.9*s.white;
% the top and the right side of the box are polygons (point lists)
hoffset = box_front(3)/5;
voffset = box_front(4)/4;
box_top = [box_front_pos(1),box_front_pos(2);
  box_front_pos(1)+box_front(3),box_front_pos(2);
  box_front_pos(1)+box_front(3)+hoffset,box_front_pos(2)-voffset;
  box_front_pos(1)+hoffset,box_front_pos(2)-voffset];
% the top is slightly darker
box_top_color = 0.85*s.white;
box_side = [box_front_pos(1)+box_front(3),box_front_pos(2);
  box_front_pos(1)+box_front(3)+hoffset,box_front_pos(2)-voffset;
  box_front_pos(1)+box_front(3)+hoffset,box_front_pos(2)+box_front(4)-voffset;
  box_front_pos(1)+box_front(3),box_front_pos(2)+box_front(4)];
% the side is darker
box_side_color = 0.8*s.white;
% make a slit on top
box_slit = [box_front_pos(1)+box_front(3)*1/4+hoffset*2/5,box_front_pos(2)-voffset*2/5;
  box_front_pos(1)+box_front(3)*3/4+hoffset*2/5,box_front_pos(2)-voffset*2/5;
  box_front_pos(1)+box_front(3)*3/4+hoffset*2/5+hoffset/5,box_front_pos(2)-voffset*3/5;
  box_front_pos(1)+box_front(3)/4+hoffset*2/5+hoffset/5,box_front_pos(2)-voffset*3/5];
% slit is dark
box_slit_color = 0.2*s.white; 
% label it as box of cards
s.box_label = 'Box of 100\nCards'; s.box_label_size = 36 / 2560 * s.winw;

s.box_front = box_front; s.box_front_pos = box_front_pos;
s.box_top = box_top; s.box_side = box_side;
s.box_front_color = box_front_color; s.box_top_color = box_top_color;
s.box_side_color = box_side_color;
s.box_slit = box_slit; s.box_slit_color = box_slit_color;
end

function s = drawBoxOfCards(s)
Screen('FillRect', s.window, s.box_front_color, s.box_front_pos);
Screen('FillPoly', s.window, s.box_top_color, s.box_top, 1);
Screen('FillPoly', s.window, s.box_side_color, s.box_side, 1);
Screen('FillPoly', s.window, s.box_slit_color, s.box_slit, 1);
oldTextSize = Screen('TextSize',s.window,s.box_label_size);
DrawFormattedText(s.window, s.box_label, 'center', 'center', ...
  s.text_color, [], [], [], [], [], s.box_front_pos);
Screen('TextSize',s.window,oldTextSize);
end

function s = cardAppear(s,im)
slit_center = mean(s.box_slit,1);
final_center = [s.winw*2/3, s.winh/2];
num_frames = floor(s.animation_time / s.ifi);
card_centers = [linspace(slit_center(1),final_center(1),num_frames);...
  linspace(slit_center(2),final_center(2),num_frames)]';
card_widths = linspace(s.box_slit(2,1)-s.box_slit(1,1),s.card_width,num_frames);
card_heights = card_widths ./ size(im,2) .* size(im,1);
card_alphas = linspace(0,1,num_frames);
% theImage = imread('./iris1.jpg');
% % Make the image into a texture
imt = Screen('MakeTexture', s.window, im);
for indf = 1:num_frames
  
  for indh = 1:numel(s.hold_on_list)
    eval(s.hold_on_list{indh});
  end
  im_rect_pos = CenterRectOnPointd([0 0 card_widths(indf) card_heights(indf)],...
    card_centers(indf,1),card_centers(indf,2));
%   Screen('PutImage', s.window, im, im_rect_pos);
  Screen('DrawTexture', s.window, imt, [], im_rect_pos,[],[], card_alphas(indf));
  Screen('Flip', s.window);
end
end

function s = cardDisappear(s,im)
slit_center = mean(s.box_slit,1);
final_center = [s.winw*2/3, s.winh/2];
num_frames = floor(s.animation_time / s.ifi);
card_centers = [linspace(slit_center(1),final_center(1),num_frames);...
  linspace(slit_center(2),final_center(2),num_frames)]';
card_widths = linspace(s.box_slit(2,1)-s.box_slit(1,1),s.card_width,num_frames);
card_heights = card_widths ./ size(im,2) .* size(im,1);
card_alphas = linspace(0,1,num_frames);
% theImage = imread('./iris1.jpg');
% % Make the image into a texture
imt = Screen('MakeTexture', s.window, im);
for indfr = 1:num_frames
  indf = num_frames+1-indfr;
  for indh = 1:numel(s.hold_on_list)
    eval(s.hold_on_list{indh});
  end
  im_rect_pos = CenterRectOnPointd([0 0 card_widths(indf) card_heights(indf)],...
    card_centers(indf,1),card_centers(indf,2));
%   Screen('PutImage', s.window, im, im_rect_pos);
  Screen('DrawTexture', s.window, imt, [], im_rect_pos,[],[], card_alphas(indf));
  Screen('Flip', s.window);
end
for indh = 1:numel(s.hold_on_list)
  eval(s.hold_on_list{indh});
end
Screen('Flip', s.window);
end

function s = makeSlider(s, y, left_text, right_text, value)

end



function testBoxOfCards()
cleanob = onCleanup(@() mycleanup);
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

s = struct;
s.hold_on_list = {};

s = makeGreyBackground(s);

s = initParams(s);

s.text_color = 0.1*s.white;

s = makeBoxOfCards(s,[s.winw/6, s.winh/3]);
s = drawBoxOfCards(s);
s.hold_on_list = [s.hold_on_list, {'s = drawBoxOfCards(s)'}];
Screen('Flip', s.window);
pause(0.2)
KbWait();
pause(0.1)
end
function testHoldOn()
cleanob = onCleanup(@() mycleanup);
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

s = struct;
s.hold_on_list = {};

s = makeGreyBackground(s);

s.text_color = 0.1*s.white;

s = makeBoxOfCards(s,[s.winw/6, s.winh/3]);
s = drawBoxOfCards(s);
s.hold_on_list = [s.hold_on_list, {'s = drawBoxOfCards(s)'}];
Screen('Flip', s.window);
pause(0.2)
KbWait();
pause(0.1)
doFadeTest()
  function doFadeTest()
    theImage = imread('./iris1.jpg');
    % Make the image into a texture
    imageTexture = Screen('MakeTexture', s.window, theImage);

    % Our will fade in and out with a sine wave function
    % See: http://en.wikipedia.org/wiki/Sine_wave
    amplitude = 0.5;
    frequency = 0.2;
    angFreq = 2 * pi * frequency;
    startPhase = 0;
    time = 0;

    % Presentation loop (press any key to exit)
    while ~KbCheck

      % Position of the square on this frame
      thisContrast = amplitude * sin(angFreq * time + startPhase) + amplitude;

      % Draw the image to the screen
      Screen('DrawTexture', s.window, imageTexture, [], [], 0, [], thisContrast);

      % draw things to hold on to the screen
      for ind = 1:numel(s.hold_on_list)
        eval(s.hold_on_list{ind});
      end
      
      % Increment the time
      time = time + s.ifi;

      % Flip to the screen
      Screen('Flip', s.window);

    end
  end
end
function testCardAppearDisappear()
cleanob = onCleanup(@() mycleanup);
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

s = struct;
s.hold_on_list = {};

s = makeGreyBackground(s);

s = initParams(s);


s = makeBoxOfCards(s,[s.winw/6, s.winh/3]);
s = drawBoxOfCards(s);
s.hold_on_list = [s.hold_on_list, {'s = drawBoxOfCards(s)'}];
Screen('Flip', s.window);

pause(0.2)
KbWait();
pause(0.1)

s = cardAppear(s, imread('./iris1.jpg'));

pause(0.2)
KbWait();
pause(0.1)

s = cardDisappear(s,imread('./iris1.jpg'));
pause(0.2)
KbWait();
pause(0.1)

end