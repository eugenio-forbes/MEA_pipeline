%%% This function updates and returns parameters used by Kilosort 
%%% based on recording information and pitch of microelectrode array
%%% Function lists signal processing and clustering parameters that
%%% can be modified.

function kilosort_configurations = n00_set_kilosort_configs(recording_directory, recording_info, MEA_pitch)

%%% Initialize structure to hold kilosort configurations
kilosort_configurations = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Modifiable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Default thresholds in kilosort3 [9 9]. Eugenio: had changed it
%%% previously to [4 2], because it wouldn't work, but now it works???
kilosort_configurations.Th = [9 9];  

%%% Number of samples to average over (annealed from first to second value)
%%% Eugenio: Default of [20 400]. Switched to [14 267] based on default
%%% sampling rate being 30 kHz.
kilosort_configurations.momentum = [14 267];

%%% Added low pass frequency in Hz  so that the data could be bandpassed 
%%% the same way as combinato extract 3000Hz or combinato detect 1000Hz
kilosort_configurations.fslow = 3000;

%%% Eugenio: Left everything below as default. %%%

%%% Time range to sort (default: all of it [0 Inf])
kilosort_configurations.trange = [0 Inf];

%%% How important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
kilosort_configurations.lam = 10;
  
%%% Splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
kilosort_configurations.AUCsplit = 0.9; 

%%% Minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
kilosort_configurations.minFR = 0.02; %Hz
 
%%% Minimum firing rate on a "good" channel (0 to skip)
kilosort_configurations.minfr_goodchannels = 0.1; %Hz

%%% Threshold crossings for pre-clustering (in PCA projection space)
kilosort_configurations.ThPre = 8;

%%% These three were main parameter changes from kilosort v2.0 to v2.5
kilosort_configurations.sig = 20;      %%% Spatial smoothness constant for registration
kilosort_configurations.fshigh = 300;  %%% High-pass frequency in Hz
kilosort_configurations.nblocks = 5;   %%% Blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.

%%% WARNING: KILOSORT RECOMMENDS NOT CHANGING ANY OF THESE OPTIONS BECAUSE OF POTENTIAL OF FATAL ERRORS !!!
%%% Options for determining PCs
kilosort_configurations.spkTh = -6;                                       %%% Spike threshold in standard deviations (-6)
kilosort_configurations.reorder = 1;                                      %%% Whether to reorder batches for drift correction. 
kilosort_configurations.nskip = 25;                                       %%% How many batches to skip for determining spike PCs
kilosort_configurations.GPU = 1;                                          %%% Has to be 1, no CPU version yet, sorry
% kilosort_configurations.Nfilt = 1024;                                   %%% Max number of clusters
kilosort_configurations.nfilt_factor = 4;                                 %%% Max number of clusters per good channel (even temporary ones)
kilosort_configurations.ntbuff = 64;                                      %%% Samples of symmetrical buffer for whitening and spike detection
kilosort_configurations.NT = 64 * 1024 + kilosort_configurations.ntbuff;  %%% Must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 

kilosort_configurations.whiteningRange = 32;                              %%% Number of channels to use for whitening each channel. 
%%% Eugenio: Default 32 channels. Could potentially change to Inf to use
%%% covariance matrix of all channels, but honestly that would make things
%%% take forever to process. Could also test increasing the number of
%%% channels by a bit to see what happens.

kilosort_configurations.nSkipCov = 25;                                    %%% Compute whitening matrix from every N-th batch
kilosort_configurations.scaleproc = 200;                                  %%% int16 scaling of whitened data
kilosort_configurations.nPCs = 3;                                         %%% How many PCs to project the spikes into
kilosort_configurations.useRAM = 0;                                       %%% Not yet available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Recording specific parameters

%%% Name of raw binary file
kilosort_configurations.fbinary = fullfile(recording_directory, 'raw_recording.bin');

%%% File name of channel map that will be used by kilosort
kilosort_configurations.chanMap = fullfile(recording_directory, 'chanMap.mat');

%%% Total number of channels in recording. Should be the sum of 'connected'.
kilosort_configurations.NchanTOT = recording_info.n_channels;

%%% File name of filtered data. Should be kept as is.
kilosort_configurations.fproc   = fullfile(recording_directory, 'temp_wh.dat');
    
%%% Sampling rate of recording in Hz
kilosort_configurations.fs = recording_info.sampling_rate;

%%% Spatial constant in um for computing residual variance of spike
%%% Eugenio: switched to actual MEA pitch from 30 um.
kilosort_configurations.sigmaMask = MEA_pitch; 

end