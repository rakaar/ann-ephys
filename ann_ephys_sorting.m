data_file_name='EPhys_data_25-Apr-2022_22_31_35-46';
 curdir=pwd;
 cd 'F:\_ANN\matalab programs\sorting\ephysdata\trial data';
eval(sprintf('load %s',data_file_name))
cd 'F:\_ANN\matalab programs\sorting';
Fs=ephysdata.Fsi;

%%% Create data matrix
[nn,kk]=size(fieldnames(ephysdata)); % f
niters=nn-1;% total iters ( how many times a stim is played )
[nsamps,nchs]=size(ephysdata.iter1);
unit_record=zeros(16,5);% what does this indicate ?



% switch PP_PARAMS.type
%     case 'AudStimEphys'
        %% General parameters
        stim_start=PP_PARAMS.protocol.stim_protocol.stim_start;% get the sim start point
        total_time=PP_PARAMS.protocol.stim_protocol.total_frames_lines;
        chnums=16;% no of channels
        %%% creating data matrix continued
        edata=zeros(niters,nsamps,chnums);
        fdata=edata;bdata=edata;
        for jj=1:niters
            eval(sprintf('edata(jj,:,:)=ephysdata.iter%i(:,1:16);',jj))
        end

        %% notch filter every channel data


        for chn=1:size(edata,3)% till the no of channels
            [fdata(:,:,chn)]=notch50filter_1ch(squeeze(edata(:,:,chn)),Fs);
        end
        % filtered data is now in fdata

        %%% baseline zeroing every channel data
        %%% fixed for 200 ms preceding stimulus
        %%% HARD CODED BASELINE DURATION
        base_p=0.2;%%% HARD CODED BASELINE DURATION                                     %%% HARD CODED BASELINE DURATION
        %%% HARD CODED BASELINE DURATION

        base_samps=fix(base_p*Fs);
        base_samps=max([1 (fix(stim_start*Fs)-base_samps)]):fix(stim_start*Fs);
        for chn=1:size(edata,3)
            [bdata(:,:,chn)]=baselinezeroing_1ch(squeeze(fdata(:,:,chn)),base_samps);
        end
        % bdata is filter+ baseline zerod data
        %% starting with analysis

        nreps=PP_PARAMS.AUD_IMG_STIM.STIMS.rep_sets;
        lps=PP_PARAMS.AUD_IMG_STIM.STIMS.lines_per_set;
        codes=PP_PARAMS.AUD_IMG_STIM.STIMS.stimcode';% take cre it should be transpose to get iters like ephysdata format 
        % got no of reps and lines per set
        SPK_Dur=3.800;%% fill with chking to sir

        %% data sorted per channel per stimuli 

        for channel_no=1:16
            fprintf('\nAalyzing Channel %i\n',channel_no);
            spike_data_channel=squeeze(bdata(:,:,channel_no));
            for stimline=1:lps
                arr_ephys_data(channel_no,stimline,:)=mean(spike_data_channel(find(codes==stimline),:));
            end
        end
randB=0;
             clc;
        %%
        for channel_no= 1:16
               fprintf('\nAalyzing Channel %i\n',channel_no);
            isonline= 0;%input(promt1);
            spike_data=squeeze(bdata(:,:,channel_no));
            for stimline=1:lps
                SPK_data(stimline,:)=mean(spike_data_channel(find(codes==stimline),:));
            end
            [no_of_units_n,no_of_units_p,ths,neg_present,pos_present,spktime_n,spktime_p,spkid_n,spkid_p,spikemat_n,spikemat_p,thres_prompt,negspcounts,posspcounts,spid_n,spid_p,spid_unit_n,spid_unit_p,spmat_n,spmat_p]=single_chnlSPKanalysismuneshwar(squeeze(bdata(:,:,channel_no)'),Fs,lps,nreps,PP_PARAMS.protocol,SPK_Dur,randB,isonline);
    threshold=thres_prompt;
                if threshold(end)==-1;
                    thresh=0;
                else
                    thresh=threshold(end-1);
                end
                curdir=pwd;
                cd sorted_data
                filestr1=strcat(data_file_name(1:end-4),'_unit_record.mat');
                if exist(filestr1)==2
                    load(filestr1)
                end
                unit_record(channel_no,:)=[thresh neg_present no_of_units_n pos_present no_of_units_p];
                unit_record_spike(channel_no)=struct('posspiketime',spktime_p,'negspiketime',spktime_n,'posspikeid',spkid_p,'negspikeid',spkid_n,'posspikemat',spikemat_p,'negspikemat',spikemat_n,'negspcounts',negspcounts,'posspcounts',posspcounts,'nspikemat_unit',spmat_n,'pspikemat_unit',spmat_p);
                filestr1=strcat(data_file_name(1:end-4),'_unit_record.mat');
                save(filestr1,'unit_record','unit_record_spike','arr_ephys_data');
                cd ..
        end 


         
