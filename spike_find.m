    function [spktime_n,spktime_p,spikemat_p,spkid_n,spikemat_n,spkid_p,n_spikes,p_spikes,offset,spike_len,rfs,stemplen,ths,thres_prompt]=spike_find(fdata_channel,Fs,lps,nreps,stimulus,SPK_Dur,randB,isonline)


% INPUT: fdata_channel=input data of single channel consisting of various
%                      iterations
%        Fs=sampling frequency of data
%        lps=lines per set i.e the no of frequency
%        nreps=no of repetion of same frequency
%        stimulus=PP_PARAMS.protocol to get stimulus data


%OUTPUT: spktime=a class containing spike instant(time) for negative as well as positive
%                threshold crossing spikes
%        spikeid=a class containing spike number at various iterations for
%                negative as well as positive spikes
%        spikemat=a matrix consisting of     positive and negative spikes of
%                 order [spike_len*spike_id]
%        spike_len=length of spikes consisting of 10 samples before
%                  threshold ccrossing and 20 after that
%        rfs=data resampled to new sampling frequency equal to 10000
%        offset=used to plot the psth at various freq decided by max value
%               in psth
%   In this function the threshold multiplier is varied and accordingly the
%   various spikes are plotted.The default threshold selected is 5 if the user
%   is not willing to change it.
%



rfs=10000; %%% frequency for resampling
[nsamples,~]=size(fdata_channel);%%%% gives no of samples and total iterations

mv=mean(fdata_channel);fdata_channel=fdata_channel-repmat(mv,nsamples,1);%%%%centering at mean i.e. subtacting mean from entire data
ffx=resample(fdata_channel,rfs,Fs);%%% resampling to 10Khz frequency
clear mv fdata_channel
%resampled_fdata_channel=ffx;
[nsamps,niters]=size(ffx);%%%% gives no of samples and total iterations after resampling
sds=std(ffx);
ths=5;last_ths=5;starting=1;transient_val=1;%% transient voltage to be eliminated starting is 1V
thres_prompt=[];
nn=1; %%% given to provide transient value just once



while last_ths~=0&last_ths~=-1
    
    if last_ths~=5 && nn==1
        promt1='give transient value(0 to exit and 1 to not give transient value)=';%%% to provide transient voltage to be removed
        transient_val= input(promt1);
        
        if transient_val~=1;
            clear spikemat_n spikemat_p spktime_n spktime_p spkid_n spkid_p ffx_v ctr_n ctr_p spikeid_n spikeid_p spike_len n_spikes p_spikes
            if transient_val==0
                ths=-1; %%% exit to channnel selection
                last_ths=-1;
            end
        else
            ths=0;
        end
        nn=0;
    end
%     if ths~=-1 && ths~=0
    prompt='ths for threshold calculation(0 to select ths and -1 to exit)=';
    ths=input(prompt);
    thres_prompt(end+1)=ths;
    
    if starting==1&ths==0
        ths=5;
    end
    if ths~=0&ths~=-1
        starting=0;
        clear spikemat_n spikemat_p spktime_n spktime_p spkid_n spkid_p ffx_v ctr_n ctr_p spikeid_n spikeid_p spike_len n_spikes p_spikes
        th=ths*sds;%%% threshold calculation
        abs_ref=0.001; %time out
        abs_ref_s=abs_ref*rfs;%%%time out in samples
        ffx_v=ffx;ctr_n=zeros(1,niters);ctr_p=zeros(1,niters);
        stemplen=0.001;stemplen_s=fix(stemplen*rfs);spbase_s=fix(stemplen_s*0.6);spikeid_n=0;spikeid_p=0;ptemplen=0.002;ptemplen_s=fix(ptemplen*rfs);spikemat_p=[];spikemat_n=[]; % changes
        %%%% creating a spike time instant for the entire iterations
        for iter_n=1:niters
            for ts=1:nsamps
                if ts-stemplen_s <= 0
                    start_t=1;
                else
                    start_t=ts-stemplen_s;
                end
                
                if ts+ptemplen_s > nsamps
                    end_t=nsamps;
                else
                    end_t=ts+ptemplen_s;
                end
                if ffx_v(ts,iter_n)>th(iter_n) && sum(ffx(start_t:end_t,iter_n) > transient_val)==0 & sum(ffx(start_t:end_t,iter_n) < -transient_val)==0 %%% second and third condition to eliminate transients
                    %positive spike occurred
                    ctr_p(iter_n)=ctr_p(iter_n)+1;
                    eval(sprintf('spktime_p.cl1.iter%i(ctr_p(iter_n))=ts/rfs;',iter_n))
                    ffx_v(ts,iter_n)=100;
                    %%%%% refractory period
                    if ts+abs_ref_s<nsamps
                        ffx_v(ts+1:ts+abs_ref_s,iter_n)=0;
                    else
                        ffx_v(ts+1:end,iter_n)=0;
                    end
                    
                elseif ffx_v(ts,iter_n)<-th(iter_n) && sum(ffx(start_t:end_t,iter_n) > transient_val)==0 & sum(ffx(start_t:end_t,iter_n) < -transient_val)==0  %%% second and third condition to eliminate transients
                    %negative spike occurred
                    ctr_n(iter_n)=ctr_n(iter_n)+1;
                    eval(sprintf('spktime_n.cl1.iter%i(ctr_n(iter_n))=ts/rfs;',iter_n))
                    ffx_v(ts,iter_n)=-100;
                    %%%%% refractory period
                    if ts+abs_ref_s<nsamps
                        ffx_v(ts+1:ts+abs_ref_s,iter_n)=0;
                    else
                        ffx_v(ts+1:end,iter_n)=0;
                    end
                end
            end
            if ctr_p(iter_n)==0
                eval(sprintf('spktime_p.cl1.iter%i=[];',iter_n))
            end
            if ctr_n(iter_n)==0
                eval(sprintf('spktime_n.cl1.iter%i=[];',iter_n))
            end
        end
        
        %%% templen is total spike time taken
        %%% ptemplen is time taken on either side of the spike time instant
        %%% crossing threshold
        %%% spkid is the class containing spike no for each iteration
        %%% spikemat is matrix of spike shapes
        %         stemplen=0.001;stemplen_s=fix(stemplen*rfs);spbase_s=fix(stemplen_s*0.6);spikeid_n=0;spikeid_p=0;ctr_p=zeros(1,niters);ctr_n=zeros(1,niters);ptemplen=0.002;ptemplen_s=fix(ptemplen*rfs);spikemat_p=[];spikemat_n=[];
        ctr_n=zeros(1,niters);ctr_p=zeros(1,niters);
        for iter_n=1:niters
            for ts=1:nsamps
                if ffx_v(ts,iter_n)==100
                    % positive spike occurred
                    ctr_p(iter_n)=ctr_p(iter_n)+1;
                    spikeid_p=spikeid_p+1;
                    eval(sprintf('spkid_p.iter%i(ctr_p(iter_n))=spikeid_p;',iter_n))
                    
                    if (ts+ptemplen_s)<=nsamps&&ts>(stemplen_s)
                        spikemat_p(:,spikeid_p)=ffx(ts-stemplen_s:ts+ptemplen_s,iter_n);
                    elseif (ts+ptemplen_s)>=nsamps %%% right side of spike exceeds the input data
                        lc=size(ffx(ts:end,iter_n),1);%%% calculation of no of samples short of
                        spikemat_p(:,spikeid_p)=[ffx(ts-stemplen_s:end,iter_n);zeros(ptemplen_s-lc+1,1)];
                    elseif ts<=(stemplen_s+1) %%%% left side of spike is short in input data
                        lc=size(ffx(1:ts,iter_n),1); %%% calculation of no of samples short of
                        spikemat_p(:,spikeid_p)=[zeros(stemplen_s-lc+1,1);ffx(1:ts+ptemplen_s,iter_n)];
                    end
                    spikemat_p(:,spikeid_p)=spikemat_p(:,spikeid_p)-mean(spikemat_p(1:spbase_s,spikeid_p));
                elseif ffx_v(ts,iter_n)==-100
                    % negative spike occurred
                    ctr_n(iter_n)=ctr_n(iter_n)+1;
                    spikeid_n=spikeid_n+1;
                    eval(sprintf('spkid_n.iter%i(ctr_n(iter_n))=spikeid_n;',iter_n))
                    
                    if (ts+ptemplen_s)<=nsamps&&ts>(stemplen_s)
                        spikemat_n(:,spikeid_n)=ffx(ts-stemplen_s:ts+ptemplen_s,iter_n);
                    elseif (ts+ptemplen_s)>=nsamps %%% right side of spike exceeds the input data
                        lc=size(ffx(ts:end,iter_n),1);%%% calculation of no of samples short of
                        spikemat_n(:,spikeid_n)=[ffx(ts-stemplen_s:end,iter_n);zeros(ptemplen_s-lc+1,1)];
                    elseif ts<=(stemplen_s+1) %%%% left side of spike is short in input data
                        lc=size(ffx(1:ts,iter_n),1); %%% calculation of no of samples short of
                        spikemat_n(:,spikeid_n)=[zeros(stemplen_s-lc+1,1);ffx(1:ts+ptemplen_s,iter_n)];
                    end
                    spikemat_n(:,spikeid_n)=spikemat_n(:,spikeid_n)-mean(spikemat_n(1:spbase_s,spikeid_n));
                end %% condition to check spike instant
            end%% ts loop
            if ctr_p(iter_n)==0
                eval(sprintf('spkid_p.iter%i=[];',iter_n))
            end
            if ctr_n(iter_n)==0
                eval(sprintf('spkid_n.iter%i=[];',iter_n))
            end
        end %% iter_n loop
        
        
        if ~isempty(spikemat_p)||~isempty(spikemat_n)%% spikemat is not empty
            
            [spike_len_n,n_spikes]=size(spikemat_n);
            [spike_len_p,p_spikes]=size(spikemat_p);
            spike_len=max(spike_len_n,spike_len_p);
            fprintf('\nGot matrix of shapes ... \n')
            
            if n_spikes==0
                fprintf('\nNo spike crossing negative threshold ... \n')
            end
            
            if p_spikes==0
                fprintf('\nNo spike crossing positive threshold ... \n')
            end
            
            find_figure('spikes');clf;
            
            %%%%%%%plot of Spikes%%%%%%%%
            
            x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
            
            if p_spikes~=0
                subplot(2,3,1),plot(x,spikemat_p,'color',[0.3 0.3 0.3])
            end
            hold on
            
            if n_spikes~=0
                subplot(2,3,1),plot(x,spikemat_n,'color',[0.6 0.6 0.6])
            end
            xlim([x(1) x(end)]);xlabel('Time(ms)');
            hold off
            clear x
            stim_start=stimulus.stim_protocol.stim_start;
            total_time=stimulus.stim_protocol.total_frames_lines;
            
            if strcmpi(stimulus.type,'rsspip')||lps==1
                offset=0;
            else
                [offset]=psth_offset(spktime_n.cl1,spktime_p.cl1,0.01,stim_start,total_time,lps,nreps,isonline,SPK_Dur);
                
            end
            if n_spikes~=0
                
                dot_raster_psth_plot(stimulus,lps,nreps,1,spktime_n,2,2,2,offset,SPK_Dur,randB,isonline);
                hold on
            end
            
            if p_spikes~=0
                dot_raster_psth_plot(stimulus,lps,nreps,1,spktime_p,2,2,1,offset,SPK_Dur,randB,isonline);
            end
            
        else
            
            fprintf('\nNo Spikes Present ... \n')
        end
        
    end % end of if
%     end
    last_ths=ths;
    if last_ths==-1
        spktime_n=[];spktime_p=[];spikemat_p=[];spkid_n=[];spikemat_n=[];spkid_p=[];n_spikes=0;p_spikes=0;offset=0;spike_len=0;stemplen=0;
    end
    
end %end of while
end

