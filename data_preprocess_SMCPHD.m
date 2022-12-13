function [GT_full,Zset,parameters,GT_per_file] = data_preprocess_SMCPHD(GT_folder,wav_folder)

    GT_files = dir(fullfile(GT_folder, '*.mat' )); 
    wav_files=dir(fullfile(wav_folder,'*.wav'));
    load(GT_files(1).name,'parameters');
    freqrange=[2000,15000];
    peak_thr=7;
    GS=0;
    Zset=[];
    
    for k=1:size(wav_files,1)
        audiofile=[wav_folder,wav_files(k).name];
        [x,fs]=audioread(audiofile);
        Zset_temp = preprocess_getZset(parameters.win_width_s,parameters.dt,fs,x,freqrange,peak_thr,GS);
        [Zset]=[Zset,Zset_temp];
    end
    
    % Make a full gt counter
    %add the individual items to the struct one by one
    
    i=1; %full GT counter
    GT_per_file=[];
    for j=1:size(GT_files,1)
            load(GT_files(j).name,'GT')
                for m=1:numel(GT)
                    GT_full(i).time=GT(m).time+60*(j-1);
                    GT_full(i).freq=GT(m).freq;
                    GT_full(i).valid=GT(m).valid;
                    i=i+1;
                end
    end

end

