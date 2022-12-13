%% /////////////////////Create lists//////////////////////////////////////
% Create a list of audio data and ground truth data files
% Make sure that ground truth data file names match the audio file names.
% Specify folder to save results to
close all;

GT_folder='./testing_process_SMCPHD/GT/Training/'
wav_folder='./testing_process_SMCPHD/audio/Training/'
GT_file =dir(GT_folder);
audio_file=dir(wav_folder);

%% ////////////// Specify PARAMETERS //////////////

%--------------Load parameters used for GT evaluation----------------------
% Note: win_width_s and dt are included in the paramaters struct
load([GT_folder,GT_file(3).name],'parameters')
training_parameters=parameters;
% ------------ Parameters for obtaining the measurements ------------------
training_parameters.freqrange=[2000,15000]; % lower and higher frequency range limits in
% Hz in which your signals occur 
training_parameters.peak_thr=training_parameters.GT_SNRthr; %threshold in dB for finding spectral peaks 
%(these become your measurements that will go to the detector)- match it to
%the threshold used to evaluate GT
training_parameters.gs=0; %Gaussian smoothing option (0 (default)= no smoothing; 1 = smoothing)
training_parameters.informed_priors=1; %0 if newborn particles are made without priors, 1 if made with knowlegde of priors
%% ////////////// TRAIN PARAMETERS from GT data//////////////

N=size(audio_file,1);
parameters_temp =cell(1,N);
for k=3:N
GTfilenameIn =[GT_folder,GT_file(k).name];
AudioFilenameIn=[audio_file(k).folder,'/',audio_file(k).name];

[parameters_temp{k}] = getParameters(GTfilenameIn,AudioFilenameIn,training_parameters);
end


detectionProbability = nan(1,N);
falseAlarmRate = nan(1,N);
birthRate = nan(1,N);
numTracks = nan(1,N);
systemNoiseDeviation = nan(1,N);
measurementDeviation = nan(1,N);
begining_freqs=[];
for k=3:N
        detectionProbability(:,k) = parameters_temp{k}.detectionProbability;
        falseAlarmRate(:,k) = parameters_temp{k}.falseAlarmRate;
        birthRate(:,k) = parameters_temp{k}.birthRate;
        numTracks(:,k) = parameters_temp{k}.numTracks;
        systemNoiseDeviation(:,k) = parameters_temp{k}.SystemNoiseDeviation;
        measurementDeviation(:,k) = parameters_temp{k}.measurementDeviation;
        begining_freqs=[begining_freqs,parameters_temp{k}.beginFreqs];
end

switch (training_parameters.informed_priors)
    case 1
        % Create a pdf based on begining frequencies
        %
        %First fit begining_freqs with different distributions using fitmethis.m
        %function:
        F= fitmethis(begining_freqs,'pdist',4); %it will plot 'pdist'=4 best 
        % distributions- make sure the one that is the best F(1) looks good

        %Compute pdf of the best distribution (F(1)- best fit):
        [y,x] = bestfitpdf(begining_freqs,F(1));
        birthpdf(1,:)=x;
        birthpdf(2,:)=y;
    case 0 
        %Uniform Distribution
        pd=makedist('Uniform','lower',training_parameters.freqrange(1),'upper',training_parameters.freqrange(2));
        x=linspace(training_parameters.freqrange(1),training_parameters.freqrange(2),101);
        birthpdf(1,:)=x
        birthpdf(2,:)=pdf(pd,x);
end

save('./testing_process_SMCPHD/Trained_parameters.mat','training_parameters',...
    'detectionProbability','falseAlarmRate','birthRate',...
    'systemNoiseDeviation','measurementDeviation','birthpdf','GT_file')


