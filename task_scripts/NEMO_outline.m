% NEMO - fMRI or behavioral task
% Latest modification: October 30, 2019
warning off

% add paths
addpath('C:\Cogent2000v1.32\Toolbox');
addpath('C:\Rosie');

% Parameters
sn = 2;
res.mode = 0;           %1 for Scan, 0 for lab, 99 for debug. res.sess should never be 4 if res.mode = 0;
res.set = 1;            % Use set 0 for odor training
res.sess = 2;           % Sess = 1-4 for set 1-4 and 1-2 for set 0.
res.start_run = 0;      % Chunk number. First run is 0.
res.end_run = 0;
res.nStim = 2;          % Number of repeats of each odor, if res.sess = 1.
res.nRep = 2;           % Number of odors presented/run (default 10) in fMRi or nunber of descriptors studied.
res.run_odors = false;  % Should use daq or not. Olfactormeter not executed if false.
displ = 0;              % 0 = window, 1= 1st screen, 2=2nd screen
delimiter = '\t';


% Derived parameters
rand('state',sum(100*clock));

chunk_path = pwd;
res.subj = sprintf('NEMO_%02d_set_%02d_sess_%02d',sn,res.set,res.sess);
res.delay = 1500;
res.sniffDur = 2000;
odor_rat = 0.9;         % Ratio of odor to clean air

if res.mode == 1
    pd = load(fullfile(chunk_path,'pd.mat'));
    res.pd = pd.pd(1:res.nRep);
end

%% ASSIGN VARIABLES TO STIMSET
% set odors
% Odor_ID, full name, [MFC1, MFC2], [Flow1, Flow2]
% [Rosie lines]
% Rosie lines 1:5 (MFC1)(6=clean air)
% Rosie lines 8:12 (MFC2)(7=clean air)
% Variables:
% example: 650, '2,3-butanedione',  [1, 7], [0.11, 0.07]
% Flow1: flow rate in MFC1
% Flow2: flow rate in M FC2
% 2,3-butanedione is in MFC1, therefore 0.11 is the flow rate for
% odor 650 and 0.07 is the flow rate for clean air

%% CONFIGURE AND START COGENT
res.textcolor = [0.2, 0.408, 0.24; 0.204, 0.224, 0.412];
config_display(displ, 3, [0.5, 0.5, 0.5], [1,1,1], 'Arial', 30, 4, 0); % cogent function. 4 indicates 4 buffers
config_keyboard;
config_mouse(10);
start_cogent;

odornames = cell(11,4);
% Instruction screen
% wait for triggerpulse

tic

for ii = res.start_run:res.end_run
    res.seed = 0+ii;                                                    % Seed for randomizing odors
    res.chunk_number = ii;                                              % Chunk ID
    chunk = file_opener_fun(res.chunk_number,chunk_path,delimiter);               % Open the CSV files
    odornames(1:end-1,:) = chunk_modifier(chunk,odor_rat,0.24);              % Backward compatibility
    odornames(end,:)= {0,'Air',[6,7],[0.12,0.12]};                      % Clean air
    res.odornames = odornames;
    
    if res.mode == 1
        settextstyle('Arial', 30);
        setforecolour(1,1,1);
        preparestring('Please wait...', 1, 0, -50);
        drawpict(1);
        waitkeydown(inf,59);                                                    % 59=enter
        clearpict(1)
        
        runairlines(res); % Pass the air
        res.select_odor = fn_select_odor_olf(res);                      % Main fMRI Task
    else
        if res.set == 0
            load(fullfile(chunk_path,sprintf('ratings%d.mat',ii)));
            res.ratings = ratings;                                      % Loaded from the ratings.mat
            res = fn_behav_training(res);                               % Odor training
        else
            settextstyle('Arial', 30);
            setforecolour(1,1,1);
            preparestring('Please wait...', 1, 0, -50);
            drawpict(1);
            waitkeydown(inf,59);                                                    % 59=enter
            clearpict(1)
            runairlines(res);
            res = fn_behav_task(res);                                   % Behavioral task
        end
    end
    
    fname = sprintf('%s_run_0%s.mat', res.subj, num2str(res.chunk_number+1));
    save(fname, 'res');
    if res.run_odors
        daq = OlfConfigDaq;                                             % Close Olfactometer, opened in fn_select_odor_olf
        OlfFlowControl(daq,0,0)
        OlfCloseLines;
    end
end
sess_time = toc;

settextstyle('Arial', 30);
setforecolour(1,1,1);
if res.mode==1
    preparestring('Please wait', 1, 0, -50);
    drawpict(1);
    waitkeydown(inf,59);                                                % 59=enter
else
    preparestring('Thank you for the participation', 1, 0, -50);
    drawpict(1);
    waitkeydown(inf,59);
end

%% STOP COGENT
stop_cogent



