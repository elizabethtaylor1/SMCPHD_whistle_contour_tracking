%Ref: https://www.mathworks.com/help/stats/clustering-using-gaussian-mixture-models.html
    
    %Input and Output data from building the RBF network
    training_data=load('InputOutput_vectors_FreqChirp.mat'); 
    %Extract Chirp Data
    X=abs(training_data.input(:,2));
    rng(3);
    k = 3; % Number of GMM components
    options = statset('MaxIter',1000);
    Sigma = {'diagonal'}; % Options for covariance matrix type
    nSigma = numel(Sigma);

    SharedCovariance = {true,false}; % Indicator for identical or nonidentical covariance matrices
    nSC = numel(SharedCovariance);
   
    for i = 1:nSigma
        for j = 1:nSC
            gmfit = fitgmdist(X,k,'CovarianceType',Sigma{i}, ...
                'SharedCovariance',SharedCovariance{j},'Options',options); % Fitted GMM
        end
    end
    save('my_gmm.mat')

