%%  T1DM_Regulator.m
%   -----------------------------------------------------------------------
%   This program simulates a set of virtual patients.The algorithms and
%   methods are based on various  mathematical models of Glucose-Insulin
%   dynamics found in literature.
%    
%   An adaptive control strategy is employed to regulate blood glucose
%   levels. This simulation platform uses an object-oriented programming
%   approach to allow for code re-usability and extra functionality.
%   Also, compartmentalization and abstraction of the processes are more
%   intuitive with this approach.
%   The controller object is also customizable to various forms to allow
%   for a testing platform.
%   -----------------------------------------------------------------------    
%   Author: Oguntuase Victor A
%   Course: Tel 799 (Research Project)
%   Thesis: Adaptive control for T1DM regulation
%   Matric No.: 202473
%   Email:  vicharde@gmail.com
%   Faculty of Technology, University of Ibadan, Nigeria.

%%
clc; clear; clear obj;       % Clearing the workspace and command window

%%
disp ('          =================== Initializing =================== ')
pause(1)
disp (' ')
disp (['      This program simulates a control test for regulating ',...
    'the BG-level' newline '      of a T1DM patient' newline])

%% Defining the simulation.

% Type of simulation ...
disp (' Please select the type of simulation to run ')
disp (' 1. for normal simulation ')
disp ([' 2. for simulation with exercise dynamics ' newline])
sim_Flag = input (' >>> ');
disp(' ')

attempt = 0;
max_NumberofTrials = 5;

while isempty(sim_Flag) && attempt < max_NumberofTrials
    attempt = attempt +1; % reducing trial attempts.
    if attempt == max_NumberofTrials
        pause(2)
        disp([newline ' You have reached maximum attempts,'...
            ' exiting now...' newline])
        pause(1)
        break
    end
    
    disp(['You have to choose a simulation type!' newline])
    pause(1)
    sim_Flag = input (' >>> ');
    disp(' ')
    
end

if sim_Flag == 1
    disp('Running simulation without exercise dynamics.')
    disp([newline 'Please wait ...'])
    pause(2)
    run virtual_Patient.m
elseif sim_Flag == 2
    disp('Running simulation with exercise dynamics.')
    disp([newline 'Please wait ...'])
    pause(2)
    run virtual_Patient2.m
end


