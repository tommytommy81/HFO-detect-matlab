% This tutorial provides the instructions to use the detector 
% described in Fedele et al. 2016, "Automatic detection of high frequency 
% oscillations during epilepsy surgery predicts seizure outcome". 
% http://www.sciencedirect.com/science/article/pii/S1388245716304394

% The first step is to create a struct "data"
% with fields
% 
% data.x_bip : matrix [Nch x samples] containing the signal in  bipolar derivation
% data.label : cell with channels name
% data.fs    : sampling frequency in Hz

% In principle any filter can be used to pass band the data. In out experience we
% had more reliable results with FIR filters. An FIR for fs = 2048 is contained in 
% FIR_2048.mat and described in details in the paper. 

% The detector can be run with the first two stages only (line 70, input.stageIII_flag  = 0), 
% the first for the baseline detection and the second for HFO validation in TF domain. 
% There is also a third stage, aimed to removed spurious events on
% the base of spatial correlation. In order to use the third stage, a
% little more effort is required by the user to implement the structure 
% Datasetup (see at line 53).

% OUTPUT of the detection
% 
% HFO_R  or HFO_FR struct with fields 
% 
% HFO_R.HFOobj:  [1x1 struct]
% HFO_R.results: [1xNch struct]
%     
% HFO_R.HFOobj contains several fields, which are commented in the code. 
% The real output is contained in HFO_R.HFOobj.EVENTS3, where each event 
% corresponds to a row, and the column are described in HFO_R.HFOobj.EVENTSlabel.
% 
% HFO_R.results is struct array 1 x Nch, with field results. 
% The same statistic contained in EVENTS3 is ordered here per channels.

% author of the code: Tommaso Fedele. 
% contact: tommaso.fedele@usz.ch


%% %%%%%%%%%%%%%%%%%%% EXAMPLE START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all, clc

%% load the data
load('patient1.mat')
data.x_bip        = EegData';
data.label        = ChanBipolar;
data.fs           = 2048;


%% create the DATASETUP struct for stage III

% this must be prepared after looking at labels.
% if you need to use it on bidimensional grid, I prepared for you various geometry
 
% the idea is to create a proximity matrix. For each contact, two matrices 
% Datasetup1(i).Dist_val: distance of each electrode from electrode i in ascending order
% Datasetup1(i).Dist_ord: id of each electrode in Dist_val according to the order in the data.x_bip and labels

%  you can find some examples for your electrode grids 
%  in the file "SETUPs.mat". For each case you have labels, Coord and Dist
%  and you need to calculate Datasetup1.

 
load SETUPs 

Coord = Setup5x4.Coord;
Dist  = Setup5x4.Dist;

for i = 1:length(Coord)    
    [Datasetup1(i).Dist_val, Datasetup1(i).Dist_ord] = sort(Dist(i,:));
end 
data.Datasetup = Datasetup1;
 
save('data', 'data')
 
 
%% Run the detector
 
load FIR_2048

load('data'); 

input.filter         = filter; % pre computed FIR
input.channel_name   = data.label; % channel labels
input.fs             = data.fs;     % sampling frequency (for TF)
input.stageIII_flag  = 0; % if 0 does not execute stage III and does not require the Datasetup struct
if input.stageIII_flag 
  input.Datasetup = data.Datasetup;
end

% %%%%%%%%%%%% ripples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.hp           = 80;       % high pass frequency
input.lp           = 250;      % low pass frequency
clear  HFOobj results
running_DetStockwell160329
HFO_R.HFOobj = HFOobj; 
HFO_R.results = results;   

% %%%%%%%%%%%% FR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.hp           = 250;      % high pass frequency
input.lp           = 500;      % low pass frequency
% input.Datasetup = Datasetup;
clear  HFOobj results
% [lm_bip_FR.HFOobj, lm_bip_FR.results]  = UMC_IIIstage_detectorFreiburg(data.x,input);
running_DetStockwell160329
HFO_FR.HFOobj = HFOobj; 
HFO_FR.results = results;     
 
save(['HFO_R_FR'],'HFO_R','HFO_FR');

 