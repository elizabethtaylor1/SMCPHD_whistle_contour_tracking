function [parameters,measurementsAssigned,measurements,measurementsDeviation,T] = getParameters(GTfilenameIn,AudioFilenameIn,training)

win_width_s=training.win_width_s;
dt=training.dt;
freqrange=training.freqrange;
peak_thr=training.peak_thr;
gs=training.gs;

% Load labeled ground truth
load(GTfilenameIn,'GT')
TracksTime = {GT.time};
TracksFreq = {GT.freq};

% Extract measurements from data
[x,fs]=audioread(AudioFilenameIn);
x=x(:,1); %select first channel
% create time vector for data based on time step:
T=dt:dt:length(x)/fs; %REPLACE THIS to be INPUT PARAM
measurements = preprocess_getZset(win_width_s,dt,fs,x,freqrange,peak_thr,gs); %In GMPHD folder
numSteps=size(measurements,2);


% Assign measurements to GT
df = 1/win_width_s; %frequency bin width
maxfreqdiff = 3*df; %maximum allowed frequency deviation from 
% annotated to measurement (needed for data association in parameter training)


numObjects = zeros(numSteps,1);
numDetections = zeros(numSteps,1);
numMissedDetections = zeros(numSteps,1);
numFalseAlarms = zeros(numSteps,1);
numBirths = zeros(numSteps,1);

measurementsSensor = cell(numSteps,1);
measurement_diff= cell(numSteps,1);

for step = 1:numSteps
    numTracks = numel(TracksTime);
    currentFreqs = zeros(0,1);
    for track = 1:numTracks
        timeDiffs = abs(TracksTime{track} - T(step)); %azimuthHistrogram.T- vector 1xNsteps- change to our time ref
        index = timeDiffs < 0.3*dt; 
        currentFreq = TracksFreq{track}(index);

        if(~isempty(currentFreq))
            currentFreqs = [currentFreqs;currentFreq];
        end

        if(index(1))
            numBirths(step,1) = numBirths(step,1) + 1;
        end
    end

    if ~isempty(measurements{step}) % if there are measurements at current time step
        currentMeasurements = measurements{step}(1,:);
        numMeasurements = size(currentMeasurements,2);

        % assign tracks to measurements
        numTracks = size(currentFreqs,1);
        weights = zeros(numMeasurements,numTracks);
        for indexM = 1:numMeasurements
            for indexT = 1:numTracks
                weights(indexM,indexT) = abs(currentMeasurements(indexM) - currentFreqs(indexT));
            end
        end
        weights(weights>maxfreqdiff) = maxfreqdiff;
        [assignmentsTmp,~] = hungarian(weights);
        assignmentsTmp(weights==maxfreqdiff) = 0;

        % extract successfully assigned measurements
        indexesM = logical(sum(assignmentsTmp,2));
        measurementsSensor{step} = currentMeasurements(indexesM);

        % Get differences between assigned measurements and tracks:
%         meas_diff=abs(repmat(currentMeasurements(:),1,numTracks) - repmat(currentFreqs',numMeasurements,1));
%         meas_diff=meas_dif(assignmentsTmp);
        meas_diff=weights(assignmentsTmp);
        if isempty(meas_diff)
            meas_diff=[]; %to ensure it's not an empty matrix
        end
        measurement_diff{step}=meas_diff;

    else % if there are no measurements at current time step
        currentMeasurements=[];
        assignmentsTmp = [];
        measurementsSensor{step} = [];
        measurement_diff{step}=[];

    end

    % save number of detections
    numObjects(step,1) = numel(currentFreqs);
    numDetections(step,1) = sum(assignmentsTmp(:));
    numMissedDetections(step,1) = numObjects(step,1) - numDetections(step,1);
    numFalseAlarms(step,1) = numel(currentMeasurements) - numDetections(step,1);

end
measurementsAssigned  = measurementsSensor;
measurementsDeviation = measurement_diff;


%% Calculate tracking parameters
    
% calculate the false alarm rate parameter as either mean or variance of the numFalseAlarms (whatever is larger)
% (note that under the Poisson assumption for false alarms used by tracking algorithms, the mean is equal to the variance
% taking the larger of the two estimated from data, makes sure that the Poison distribution is conservative and not overconfident)
parameters.falseAlarmRate = max(var(numFalseAlarms),mean(numFalseAlarms));

% calculate detection probability by averaging over all time steps with a least one object present
detectionProbabilitySteps = numDetections./numObjects;
detectionProbabilitySteps = detectionProbabilitySteps(~isnan(detectionProbabilitySteps));
parameters.detectionProbability = mean(detectionProbabilitySteps);

parameters.numTracks = numel(TracksFreq);
parameters.birthRate = max(var(numBirths),mean(numBirths));

% Compute measurement deviation
parameters.measurementDeviation=mean(cellfun(@mean,measurement_diff),'omitnan');


% calculate maximum acceleration- for estimating system noise deviation
%dt = (currentTracksTimeSensor{1}(2)-TracksTime{1}(1)); %time step in seconds. Assumes time steps are uniform throughout
accelerations = cellfun(@(x) diff(diff(x)'./dt)./dt,TracksFreq,'UniformOutput',false);
maxAccelerationTmp = cellfun(@(x) max(abs(x)),accelerations,'UniformOutput',false); %maximum acceleration per track
maxAcceleration = vertcat(maxAccelerationTmp{:});


parameters.SystemNoiseDeviation = median(maxAcceleration); % Take median of the maximum accelerations to be the typicall max acceleration
% Take SystemNoiseDeviation to be 0.5*a_M < SystemNoiseDeviation < a_M

parameters.beginFreqs = cellfun(@(x) x(1),TracksFreq);


end