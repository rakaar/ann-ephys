function muneshwar(data_file_path,data_file_name)
%and all spike time, spike id,spike matrix, spike count will be
%stored as structure as unit_record_spike file in same unit_analysis folder
curr_dir=pwd;
cd(data_file_path);
%eval(sprintf('cd %s',data_file_path))
eval(sprintf('load %s',data_file_name(1:end-4)))
cd(curr_dir);
%eval(sprintf('cd %s',curr_dir))
Fs=ephysdata.Fsi;

%%% Create data matrix
[nn,kk]=size(fieldnames(ephysdata)); % fields = fieldnames(S) returns the field names of the structure array S in a cell array.
niters=nn-1;% bcz last field s of fsi
[nsamps,nchs]=size(ephysdata.iter1);
unit_record=zeros(16,5);

switch PP_PARAMS.type
    case 'AudStimEphys'
        %%% General parameters
        stim_start=PP_PARAMS.protocol.stim_protocol.stim_start;% get the sim start point
        total_time=PP_PARAMS.protocol.stim_protocol.total_frames_lines;

        if PP_PARAMS.protocol.mult_channel==1 %%% assuming 16 channel data
            chnums=16;% no of channels
            %%% creating data matrix continued
            edata=zeros(niters,nsamps,chnums);fdata=edata;bdata=edata;
            for jj=1:niters
                eval(sprintf('edata(jj,:,:)=ephysdata.iter%i(:,1:16);',jj))
            end
            clear ephysdata % y this?
            %%% notch filter every channel data
            for chn=1:size(edata,3)
                [fdata(:,:,chn)]=notch50filter_1ch(squeeze(edata(:,:,chn)),Fs);
            end

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

            %%% SPIKE ANALYSIS FOR SINGLE CHANNEL
            nreps=PP_PARAMS.AUD_IMG_STIM.STIMS.rep_sets;
            lps=PP_PARAMS.AUD_IMG_STIM.STIMS.lines_per_set;
            run_channels=1;
            %%%%%%%%%%%

            %%%%%%%%%%%%

            for channel_no=1:16
                fprintf('\nAalyzing Channel %i\n',channel_no);

                %%% Only for display purposes here

                %%%% code for random position of stimulus B in SSA %%%%

                if strcmpi(PP_PARAMS.protocol.type,'SSA')
                    randB=PP_PARAMS.protocol.stim_protocol.randB_check;
                else
                    randB=0;
                end
                %                 promt1=e'Is it Online analysis or Offline analysis( input 1 for online and 0 for offline)=';
                isonline= 0;%input(promt1);

                SPK_Dur=.800;  %%% 800 msec only for random SSA
                if isonline && randB==1

                    stimB_loc = PP_PARAMS.protocol.stim_protocol.rand_stimB_loc;
                    stimB_loc_1=find(stimB_loc==1);
                    spike_data=squeeze(bdata(:,:,channel_no));


                    if ~isempty(stimB_loc_1)

                        for i=length(stimB_loc_1):-1:1
                            eval(sprintf('spike_data(stimB_loc_1(%i):lps:end,:)=[];',i))
                            eval(sprintf('stimB_loc(stimB_loc_1(%i))=[];',i))
                        end
                    end

                    SPK_stimB=zeros(length(stimB_loc),SPK_Dur*Fs);
                    SPK_stimA=zeros(length(stimB_loc),SPK_Dur*Fs);

                    for ii=1:length(stimB_loc)
                        eval(sprintf('stimB_pos=stimB_loc(%i);',ii))
                        stimB_start=stim_start+(PP_PARAMS.protocol.stim_protocol.gap_dur+PP_PARAMS.protocol.stim_protocol.stim_Adur)*(stimB_pos-1);
                        if (stimB_start+SPK_Dur)*Fs > size(spike_data,2)
                            eval(sprintf('SPK_stimB(%i,:)=spike_data(%i,int64(stimB_start*Fs)+1:end);',ii,ii))
                        else
                            eval(sprintf('SPK_stimB(%i,:,:)=spike_data(%i,int64(stimB_start*Fs)+1:int64((stimB_start+SPK_Dur)*Fs));',ii,ii))
                        end
                    end
                    SPK_stimA=spike_data(:,int64(stim_start*Fs)+1:int64((stim_start+SPK_Dur)*Fs),:);
                    SPK_data=vertcat(SPK_stimA,SPK_stimB)';
                    [no_of_units_n,no_of_units_p,ths,neg_present,pos_present,spktime_n,spktime_p,spkid_n,spkid_p,spikemat_n,spikemat_p,thres_prompt,negspcounts,posspcounts,spid_n,spid_p,spid_unit_n,spid_unit_p,spmat_n,spmat_p]=single_chnlSPKanalysismuneshwar(SPK_data,Fs,length(stimB_loc),nreps*2,PP_PARAMS.protocol,SPK_Dur,randB,isonline);
                else
                    [no_of_units_n,no_of_units_p,ths,neg_present,pos_present,spktime_n,spktime_p,spkid_n,spkid_p,spikemat_n,spikemat_p,thres_prompt,negspcounts,posspcounts,spid_n,spid_p,spid_unit_n,spid_unit_p,spmat_n,spmat_p]=single_chnlSPKanalysismuneshwar(squeeze(bdata(:,:,channel_no)'),Fs,lps,nreps,PP_PARAMS.protocol,SPK_Dur,randB,isonline);
                end
                threshold=thres_prompt;
                if threshold(end)==-1;
                    thresh=0;
                else
                    thresh=threshold(end-1);
                end
                %%% for error cases makinf structure
                cd sorted_data
                filestr1=strcat(data_file_name(1:end-4),'_unit_record.mat');
                if exist(filestr1)==2
                    load(filestr1)
                end
                unit_record(channel_no,:)=[thresh neg_present no_of_units_n pos_present no_of_units_p];
                unit_record_spike(channel_no)=struct('posspiketime',spktime_p,'negspiketime',spktime_n,'posspikeid',spkid_p,'negspikeid',spkid_n,'posspikemat',spikemat_p,'negspikemat',spikemat_n,'negspcounts',negspcounts,'posspcounts',posspcounts,'nspikemat_unit',spmat_n,'pspikemat_unit',spmat_p);
                %%% change for channel_no

                %
                filestr1=strcat(data_file_name(1:end-4),'_unit_record.mat');
                %                 eval(sprintf('load %s;',filestr1))
                %                 unit_record_spike(channel_no).nspikemat_unit=spmat_n;
                %                 unit_record_spike(channel_no).pspikemat_unit=spmat_p;
                %                 unit_record_spike(channel_no).nspikemat_unit.cl1=negspikemat;
                %                 unit_record_spike(channel_no).pspikemat_unit=0;
                save(filestr1,'unit_record','unit_record_spike')
                cd ..
            end
        end

end
% cd E:\Analysis\check
% filestr1=strcat(data_file_name(1:end-4),'_unit_record.mat');
% save(filestr1,'unit_record_spike1')
% cd E:\Analysis\Unit_Analysis_Results
% filestr1=strcat(data_file_name(1:end-4),'_unit_record.mat');
% eval(sprintf('load %s;',filestr1))
% unit_record_spike(channel_no).nspikemat_unit=spmat_n;
% unit_record_spike(channel_no).pspikemat_unit=spmat_p;
% save(filestr1,'unit_record','unit_record_spike','results_spontrate','results_offset','results_response')
cd ..
end

