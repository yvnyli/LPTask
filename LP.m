%%%%%%% this part of the script gets executed

% for a mac, uncomment the line below to ignore the synchronization problem
% Screen('Preference', 'SkipSyncTests', 1);   

% run the experiment
% results = LPTask();


testTrySlider();






%%%%%%% everything below are functions, some are utility functions, some
%%%%%%% are procedural blocks organized in functions

function LPTask()
% the main function in charge of running the entire experiment

% Clear the workspace
close all;
clearvars;
sca;
cleanob = onCleanup(@() mycleanup); % see MATLAB's onCleanup
defaultSetup();

% every about the experiment is communicated between functions through the
% struct as inputs and outputs
s = struct;

s = makeGreyBackground(s);

s = initParams(s);

% s = practice(s);




s = makeBoxOfCards(s,[s.winw/6, s.winh/3]);

% s = drawBoxOfCards(s);
% s = addHoldOn(s,'s = drawBoxOfCards(s);');



end





%%%%% procedural functions

function s = practice(s)
% participant receives instructions, learns about keyboard/mouse controls,
% and do practice round(s)
end

function s = trySlider(s)
% participant learns how to use the slider here
% this can be a subroutine in practice()

% slider for demonstration
[s,sind] = initSlider(s,s.winh*2/3,'0','100',50,...
  'Please drag the slider to 100.',1,@returnDown);

% main instruction
textbox = CenterRectOnPointd([0 0 100 100],s.xc,s.winh/3);
[s,hind1] = addTextHoldOn(s,[],s.big_text_size,s.window,...
  'When the slider is red, it is active and you can change its value with the mouse.',...
  'center','center',s.text_color,[],[],[],[],[],textbox);

% move-on instruction
[s,hind2] = addTextHoldOn(s,[],s.bottom_text_size,s.window,...
  'Press ENTER to confirm your response.',...
  'center','center',s.text_color,[],[],[],[],[],s.bottom_text_box);

[s,hind3] = addSliderHoldOn(s,sind,[]);
s = activateSlider(s,sind);

s = deleteHoldOn(s,hind2); s = deleteHoldOn(s,hind1); hind3 = hind3-2;

DrawFormattedTextPlus(s.big_text_size,s.window,...
  'Once you confirm your response, the slider can''t be changed.',...
  'center','center',s.text_color,[],[],[],[],[],textbox);
DrawFormattedTextPlus(s.bottom_text_size,s.window,...
  'Press SPACEBAR to continue.',...
  'center','center',s.text_color,[],[],[],[],[],s.bottom_text_box);

drawHoldOn(s);
Screen('Flip',s.window);

s = deleteHoldOn(s,hind3); 
while 1
  if spaceDown()
    break;
  end
end
s = deleteSlider(s,sind);
end

function s = runTrial(s,isPractice)
% this function runs a trial (practice or for real) including saving the
% data

if isPractice
  % sample 2 and then 2
else
  % sample according to task design
end
end












%%%%%% utility functions

function mycleanup()
% when the program crashes or when it finishes
sca; % get rid of the screen
ListenChar(0); % give keyboard control back to MATLAB
end

function defaultSetup()
% all the setup things
PsychDefaultSetup(2); % what i'm supposed to do
KbName('UnifyKeyNames'); % because why not
ListenChar(2); % prevents keypresses to be entered into MATLAB
end

function s = makeGreyBackground(s)
% sets up the "screen" which is a grey background and gets all the
% parameters about the moniter/window/screen

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
% Maximum priority level for better-controlled timing
topPriorityLevel = MaxPriority(s.window);
Priority(topPriorityLevel);

% Query the frame duration
s.ifi = Screen('GetFlipInterval', s.window);

% Get the centre coordinate and dimensions of the window
[s.xc, s.yc] = RectCenter(windowRect);
s.winw = windowRect(3);
s.winh = windowRect(4); 

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', s.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
end

function s = initParams(s)
% this is where all the parameters should be defined for better
% organization

s.big_text_size = 36  / 2560 * s.winw;
s.bottom_text_size = 36 / 2560 * s.winw;

s.card_width = s.winw/5; % size of stimulus
s.animation_time = 0.5; % duration of in/out animation

s.text_color = 0.1*s.white;
s.text_color_active = 0.9*s.white; % for activated slider

s.slider_color = 0.2*s.white; % regular grey slider tick and scale
s.slider_color_active = [0.9*s.white s.black s.black]; % red slider tick

s.slider_len = s.winw/4; % length of the slider scale
s.slider_thickness = s.slider_len / 512; % thickness of slider scale
s.slider_w = s.slider_len / 100; % width of slider tick
s.slider_h = s.slider_w * 5; % height of slider tick
s.slider_text_size = 36 / 2560 * s.winw; % size of slider labels

s.text_dummy_box = [0 0 100 100]; % for using DrawFormattedText
% placement of the move-on instructions (e.g. Press ENTER)
s.bottom_text_box = CenterRectOnPointd([0 0 100 100],s.xc,s.winh*7/8);


% hold on list is a list of things to hold onto the screen. it's a list of
% commands (i.e. texts, hence the cell array) that draw those things
% this allows the function that creates things to be the only one
% responsible for those things, and its subroutines can be agnostic to
% what's on the screen
% see addHoldOn() for more 
s.hold_on_list = {}; 


% sliders parameters are kept as structs and in this array for easy
% re-drawing
s.sliders = [];


end

function s = makeBoxOfCards(s, front_center)
% make a box of cards centered on a point and save components to s for
% re-drawing
% front_center: the center coord of the front size of the box
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
% save thing to the struct
s.box_front = box_front; s.box_front_pos = box_front_pos;
s.box_top = box_top; s.box_side = box_side;
s.box_front_color = box_front_color; s.box_top_color = box_top_color;
s.box_side_color = box_side_color;
s.box_slit = box_slit; s.box_slit_color = box_slit_color;
end

function s = drawBoxOfCards(s)
% assuming the box of cards has been made, this function draw it to the
% screen
Screen('FillRect', s.window, s.box_front_color, s.box_front_pos);
Screen('FillPoly', s.window, s.box_top_color, s.box_top, 1);
Screen('FillPoly', s.window, s.box_side_color, s.box_side, 1);
Screen('FillPoly', s.window, s.box_slit_color, s.box_slit, 1);
DrawFormattedTextPlus(s.box_label_size, s.window, s.box_label, ...
  'center', 'center', s.text_color, [], [], [], [], [], s.box_front_pos);
end

function s = cardAppear(s,im)
% do the animation of a card appearing
% im: the image returned by MATLAB's imread
% the card comes out of the center of the slit of the box
slit_center = mean(s.box_slit,1);
% this is where its center will be after the animation
final_center = [s.winw*2/3, s.winh/2];
% calculate the position, dimension, and transparency of the card in each frame
num_frames = floor(s.animation_time / s.ifi);
card_centers = [linspace(slit_center(1),final_center(1),num_frames);...
  linspace(slit_center(2),final_center(2),num_frames)]';
card_widths = linspace(s.box_slit(2,1)-s.box_slit(1,1),s.card_width,num_frames);
card_heights = card_widths ./ size(im,2) .* size(im,1);
card_alphas = linspace(0,1,num_frames);

% im = imread('./iris1.jpg'); % assuming this has happened
% Make the image into a texture
imt = Screen('MakeTexture', s.window, im);
for indf = 1:num_frames  
  drawHoldOn(s);
  im_rect_pos = CenterRectOnPointd([0 0 card_widths(indf) card_heights(indf)],...
    card_centers(indf,1),card_centers(indf,2));
%   Screen('PutImage', s.window, im, im_rect_pos); % this is too slow
  Screen('DrawTexture', s.window, imt, [], im_rect_pos,[],[], card_alphas(indf));
  Screen('Flip', s.window);
end
end

function s = cardDisappear(s,im)
% animation of a card disappearing
% im: the image returned by MATLAB's imread; caller responsible for using
% the right image and calling appear before disappear
% same exact thing as cardAppear() except the index of the loop is in
% reverse
slit_center = mean(s.box_slit,1);
final_center = [s.winw*2/3, s.winh/2];
num_frames = floor(s.animation_time / s.ifi);
card_centers = [linspace(slit_center(1),final_center(1),num_frames);...
  linspace(slit_center(2),final_center(2),num_frames)]';
card_widths = linspace(s.box_slit(2,1)-s.box_slit(1,1),s.card_width,num_frames);
card_heights = card_widths ./ size(im,2) .* size(im,1);
card_alphas = linspace(0,1,num_frames);
imt = Screen('MakeTexture', s.window, im);
for indfr = 1:num_frames
  indf = num_frames+1-indfr; % reverse the index
  drawHoldOn(s);
  im_rect_pos = CenterRectOnPointd([0 0 card_widths(indf) card_heights(indf)],...
    card_centers(indf,1),card_centers(indf,2));
  Screen('DrawTexture', s.window, imt, [], im_rect_pos,[],[], card_alphas(indf));
  Screen('Flip', s.window);
end
drawHoldOn(s);
Screen('Flip', s.window);
end

function [s,sind] = initSlider(s,y,left_text,right_text,value,...
  prompt_text,numbered,stoppingFunc)
% initialize a slider with the given parameters
% y: y position of the slider
% left_text: text on the left end of the slider
% right_text: text on the right end of the slider
% value: initial value of the slider
% prompt_text: text of the prompt of the slider (displayed above slider)
% numbered: 0 or 1; whether the value of the slider is displayed
% stoppingFunc: the slider sequence will stop when the stoppingFunc returns
%   true. this is so that the creator of the slider can decide what to use
%   without altering the function that runs the sequence
% sind: slider index in the array s.sliders
sl.y = y; sl.left_text = left_text; sl.right_text = right_text;
sl.value = value; sl.prompt_text = prompt_text; sl.numbered = numbered;
sl.stoppingFunc = stoppingFunc;
s.sliders = [s.sliders, sl];
sind = numel(s.sliders);
end

function s = deleteSlider(s,sind)
% delete the sind'th slider
s.sliders(sind) = [];
end

function s = makeAndDrawSlider(s, sind)
% make components of the sind'th slider and draw it (inactivated) to the screen
slider_bar = [0 0 s.slider_len, s.slider_thickness];
slider_bar_pos = CenterRectOnPointd(slider_bar, s.xc, s.sliders(sind).y);
Screen('FillRect', s.window, s.slider_color, slider_bar_pos);
text_dummy_box_left = CenterRectOnPointd(s.text_dummy_box,...
  slider_bar_pos(1),slider_bar_pos(2)+s.text_dummy_box(4)/2);
text_dummy_box_right = CenterRectOnPointd(s.text_dummy_box,...
  slider_bar_pos(3),slider_bar_pos(2)+s.text_dummy_box(4)/2);
text_dummy_box_mid = CenterRectOnPointd(s.text_dummy_box,...
  s.xc,slider_bar_pos(2)-s.text_dummy_box(4));
DrawFormattedTextPlus(s.slider_text_size,s.window, s.sliders(sind).left_text,...
  'center', 'center', s.text_color, [], [], [], [], [], text_dummy_box_left);
DrawFormattedTextPlus(s.slider_text_size,s.window, s.sliders(sind).right_text,...
  'center', 'center', s.text_color, [], [], [], [], [], text_dummy_box_right);
DrawFormattedTextPlus(s.slider_text_size,s.window, s.sliders(sind).prompt_text,...
  'center', 'center', s.text_color, [], [], [], [], [], text_dummy_box_mid);
if s.sliders(sind).numbered
  text_dummy_box_mid(2) = slider_bar_pos(2)-s.text_dummy_box(4)/2;
  DrawFormattedTextPlus(s.slider_text_size,s.window, num2str(s.sliders(sind).value),...
    'center', 'center', s.text_color, [], [], [], [], [], text_dummy_box_mid);
end
tick_pos = CenterRectOnPointd([0 0 s.slider_w, s.slider_h], ...
  slider_bar_pos(1) + s.slider_len * s.sliders(sind).value / 100, s.sliders(sind).y);
% s.slider_color is the inactive color
Screen('FillRect', s.window, s.slider_color, tick_pos);
end

function s = activateSlider(s, sind)
% loop for handling slider controls including
%   1. click on the slider scale: slider jumps to position
%   2. hold and drag
% sind: index of slider in s.sliders
slider_bar = [0 0 s.slider_len, s.slider_thickness];
slider_bar_pos = CenterRectOnPointd(slider_bar, s.xc, s.sliders(sind).y);
% area of the slider bar but of the same height as slider tick
slider_area = CenterRectOnPointd([0 0 s.slider_len s.slider_h], s.xc, s.sliders(sind).y);
vbl = Screen('Flip', s.window);
holding = false;
while 1
  % loop control until stopping criterion is met (e.g. ENTER is pressed)
  % this is set in initSlider
  if s.sliders(sind).stoppingFunc()
    break;
  end

  [mx,my,buttons] = GetMouse(s.window);
  inside = IsInRect(mx,my,slider_area); % whether mouse is inside slider area
  holding = holding && any(buttons); % if buttons are not pressed, end holding status
  if (inside && any(buttons)) || holding 
    % if mouse is inside and button is pressed, or the button is held down
    % (regardless of where the mouse is), change slider value according to
    % mouse position
    s.sliders(sind).value = ...
      min(max(round((mx - slider_bar_pos(1))/s.slider_len*100),0),100);
   end
  if ~holding
    % to start the "holding" status, mouse must be inside the scale area
    holding = any(buttons) && inside;
  end
  % now redraw everything, note that the slider (inactive version) itself
  % is in the hold on list
  drawHoldOn(s);
  % the only things needed to be drawn are the red tick and white texts
  % (covering grey tick and grey texts)
  tick_pos = CenterRectOnPointd([0 0 s.slider_w, s.slider_h], ...
    slider_bar_pos(1) + s.slider_len * s.sliders(sind).value / 100, s.sliders(sind).y);
  Screen('FillRect', s.window, s.slider_color_active, tick_pos);
  text_dummy_box_mid = CenterRectOnPointd(s.text_dummy_box,...
    s.xc,slider_bar_pos(2)-s.text_dummy_box(4));
  oldTextSize = Screen('TextSize',s.window,s.slider_text_size);
  DrawFormattedText(s.window, s.sliders(sind).prompt_text, 'center', 'center', ...
    s.text_color_active, [], [], [], [], [], text_dummy_box_mid);
  if s.sliders(sind).numbered
    text_dummy_box_mid(2) = slider_bar_pos(2)-s.text_dummy_box(4)/2;
    DrawFormattedText(s.window, num2str(s.sliders(sind).value), 'center', 'center', ...
      s.text_color_active, [], [], [], [], [], text_dummy_box_mid);
  end
  Screen('TextSize',s.window,oldTextSize);

  vbl  = Screen('Flip', s.window, vbl + (1 - 0.5) * s.ifi);
end
% after loop ends, redraw so that the slider appears inactive
drawHoldOn(s);
Screen('Flip', s.window);
end

%%%%%% a couple of candidate stoppingFuncs

function tf = returnDown()
[~,~,keycode,~] = KbCheck();
tf =  keycode(KbName('Return'));
end

function tf = spaceDown()
[~,~,keycode,~] = KbCheck();
tf =  keycode(KbName('space'));
end

function tf = downArrowDown()
[~,~,keycode,~] = KbCheck();
tf =  keycode(KbName('downarrow'));
end

function tf = upArrowDown()
[~,~,keycode,~] = KbCheck();
tf =  keycode(KbName('uparrow'));
end

function [s,hind] = addHoldOn(s,command,afterwhich)
% adds a command to the hold on list for re-drawing a thing everytime the
% screen flips, see companion function drawHoldOn()
% command: the command to be run to draw the thing
% afterwhich: the item will be inserted after the afterwhich'th command,
%   if this is left blank (as []), it is defaulted to the last item
% hind: the index of the new item in the hold on list

if isempty(afterwhich) || afterwhich>=numel(s.hold_on_list)
  % insert to the end
  s.hold_on_list = [s.hold_on_list, {command}];
  hind = numel(s.hold_on_list);
else
  % insert after afterwhich
  s.hold_on_list = [s.hold_on_list(1:afterwhich),{command},s.hold_on_list((afterwhich+1):end)];
  hind = afterwhich + 1;
end
end

function [s,hind] = addTextHoldOn(s,afterwhich,textsize,...
  ~,tstring,sx,sy,color,wrapat,fH,fV,vS,rtl,winRect)
% adds a formatted text to hold on, this is a wrapper on DrawFormattedText
% afterwhich: same as in addHoldOn()
% textsize: size of the text! can't believe its not part of
%   DrawFormattedText already
% ~,...,winRect are the same as parameters to DrawFormattedText()
% (~ is win but we always use s.window)

% the trick is to convert everything into text
if isempty(sx) % if it is [], make it '[]'
  sxtext = '[]';
elseif isnumeric(sx) % if it is a number, make it string
  sxtext = num2str(sx);
else % else it must be text, add '' around it
  sxtext = ['''',sx,''''];
end
if isempty(sy) % same as sx for sy
  sytext = '[]';
elseif isnumeric(sy)
  sytext = num2str(sy);
else
  sytext = ['''',sy,''''];
end
if isempty(fH) % flip horizontal: if it is [], put default value 0
  fH = 0;
end
if isempty(fV) % flip vertical: same as fH
  fV = 0;
end
if isempty(vS) % vSpacing: if it is [], put default value 1
  vS = 1;
end
if isempty(rtl) % righttoleft: same as fH
  rtl = 0;
end
% make sure these vectors/scalars are horizontal
color = color(:)'; winRect = winRect(:)'; 
% DrawFormattedText
dft = ...
  sprintf('DrawFormattedText(s.window,''%s'',%s,%s,[%s],[%s],%d,%d,%d,%d,[%s]);',...
  tstring,sxtext,sytext,num2str(color),num2str(wrapat),fH,fV,vS,rtl,num2str(winRect));
% change text size
change_size = ...
  sprintf('oldTextSize = Screen(''TextSize'',s.window,%.2f);',textsize);
% change text size back
change_back = 'Screen(''TextSize'',s.window,oldTextSize);';
% compile the above 3 commands and add it to hold on
[s,hind] = addHoldOn(s,[change_size,dft,change_back],afterwhich);
end

function DrawFormattedTextPlus(textsize,...
  win,tstring,sx,sy,color,wrapat,fH,fV,vS,rtl,winRect)
% basically DrawFormattedText but with the option of changing text size
oldTextSize = Screen('TextSize',win,textsize);
DrawFormattedText(win,tstring,sx,sy,color,wrapat,fH,fV,vS,rtl,winRect);
Screen('TextSize',win,oldTextSize);
end

function [s,hind] = addSliderHoldOn(s,sind,afterwhich)
% add slider to hold on list
% sind: the index of the slider in the s.sliders to draw
% afterwhich: passed into addHoldOn
% hind: the output index of addHoldOn
[s,hind] = addHoldOn(s,['makeAndDrawSlider(s,',num2str(sind),');'],afterwhich);
end

function drawHoldOn(s)
% go through the command in the hold on list and draw all of them
for indh = 1:numel(s.hold_on_list)
  eval(s.hold_on_list{indh});
end
end

function s = deleteHoldOn(s,hind)
% delete the hind'th hold on thing from the list
s.hold_on_list(hind) = [];
end






%%%% as I make features, I test them with the functions below, and they are
%%%% not for production use

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
function testmakeAndDrawSlider()
cleanob = onCleanup(@() mycleanup);
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

s = struct;

s = makeGreyBackground(s);

s = initParams(s);


s = makeBoxOfCards(s,[s.winw/6, s.winh/3]);
s = drawBoxOfCards(s);
s.hold_on_list = [s.hold_on_list, {'s = drawBoxOfCards(s);'}];
Screen('Flip', s.window);

pause(1)

drawHoldOn(s);
[s,sind] = initSlider(s,1000,'not at all','completely',42,'How confident are you?',true);
s = makeAndDrawSlider(s,sind);
s.hold_on_list = [s.hold_on_list, {['s = makeAndDrawSlider(s,', num2str(sind),');']}];
Screen('Flip', s.window);

pause(1)

s = cardAppear(s, imread('./iris1.jpg'));

pause(1)

s = cardDisappear(s,imread('./iris1.jpg'));
pause(0.2)
KbWait();
pause(0.1)

end
function testActivateSlider()
cleanob = onCleanup(@() mycleanup);
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
KbName('UnifyKeyNames');

s = struct;

s = makeGreyBackground(s);

s = initParams(s);


s = makeBoxOfCards(s,[s.winw/6, s.winh/3]);
s = drawBoxOfCards(s);
s.hold_on_list = [s.hold_on_list, {'s = drawBoxOfCards(s);'}];
Screen('Flip', s.window);

pause(1)

drawHoldOn(s);
[s,sind] = initSlider(s,1000,'not at all','completely',42,'How confident are you?',true);
s = makeAndDrawSlider(s,sind);
s.hold_on_list = [s.hold_on_list, {['s = makeAndDrawSlider(s,', num2str(sind),');']}];
Screen('Flip', s.window);

pause(0.1)
KbWait();

s = activateSlider(s,sind);

pause(0.1);
KbWait();
pause(0.1);
end
function testTrySlider()
cleanob = onCleanup(@() mycleanup);
defaultSetup();
s = struct;

s = makeGreyBackground(s);

s = initParams(s);

s = trySlider(s);
pause(0.1);KbWait();
end


%%%%%%% end of file