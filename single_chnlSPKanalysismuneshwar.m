function [k_mean_n,k_mean_p,ths,neg_spikes,pos_spikes,spktime_n,spktime_p,spkid_n,spkid_p,spikemat_n,spikemat_p,thres_prompt,spcounts_n,spcounts_p,spid_n,spid_p,spid_unit_n,spid_unit_p,spmat_n,spmat_p]=  single_chnlSPKanalysismuneshwar(fdata_channel,Fs,lps,nreps,stimulus,SPK_Dur,randB,isonline)

k_mean_p=0;pos_spikes=0;
k_mean_n=0;neg_spikes=0;
spcounts_n=0;spcounts_p=0;
spid_n=0;spid_p=0;
spid_unit_p=0;spid_unit_n=0;
spmat_n=0;spmat_p=0;
[spktime_n,spktime_p,spikemat_p,spkid_n,spikemat_n,spkid_p,n_spikes,p_spikes,offset,spike_len,rfs,stemplen,ths,thres_prompt]=spike_find(fdata_channel,Fs,lps,nreps,stimulus,SPK_Dur,randB,isonline);

if ths~=-1
    if ~isempty(spikemat_n)
        if p_spikes~=0
            w = input('Select 1 to include positive threshold crossing or select 2 for only negative threshold crossing (1,2)?','s');
            if strcmpi(w,'1')
                neg_spikes=1;
                pos_spikes=1;
                
                remove_spikes_n = input('Do you want to do remove  negative crossing spikes before clustering (y/n)?','s');
                
                if strcmpi(remove_spikes_n,'y')
                    idx=zeros(1,size(spikemat_n,2)); spikemat=[]; spkid=[]; spktime=[]; n_spikes =[];
                    spikemat = spikemat_n;
                    spkid = spkid_n;
                    spktime = spktime_n.cl1;
                    spikemat_n=[]; spkid_n=[]; spktime_n.cl1=[];
                    
                   [spikemat_n,spkid_n,spktime_n.cl1,~]= spike_removal_gui(spikemat,spkid,spktime,idx);
                   n_spikes = size(spikemat_n,2);
                end
                
                remove_spikes_p = input('Do you want to do remove positive crossing spikes before clustering (y/n)?','s');
                
                if strcmpi(remove_spikes_p,'y')
                    idx=zeros(1,size(spikemat_p,2)); spikemat=[]; spkid=[]; spktime=[]; p_spikes =[];
                    spikemat = spikemat_p;
                    spkid = spkid_p;
                    spktime = spktime_p.cl1;
                    spikemat_p=[]; spkid_p=[]; spktime_p.cl1=[];
                    
                   [spikemat_p,spkid_p,spktime_p.cl1,~]= spike_removal_gui(spikemat,spkid,spktime,idx);
                   p_spikes = size(spikemat_p,2);
                end
                
                find_figure('Selection of K for k-mean clustering');clf;
                find_figure('Selection of nc');clf;
                
                
                %%%% for negative spikes%%%%%%%
                
                [nc_n,k_mean_n,k_centres_n]=nc_kmean(spikemat_n,n_spikes,p_spikes,2,1,0,0,0);
                [spcounts_n,spid_n,spid_unit_n,spmat_n,centers_n]=class_sorting(spikemat_n,spktime_n.cl1,spkid_n,nc_n,k_mean_n,2,0,0,k_centres_n);
                v = input('Do you want to do manual clustering (y/n)?','s');
                
                if strcmpi(v,'y')
                    spcounts_n=[]; centers_n=[]; spmat_n=[];
                    [spcounts_n,centers_n,spmat_n]= manual_cluster(spikemat_n,spkid_n,spktime_n.cl1);
                end
                find_figure('spikes');clf;
                
                
                x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
                subplot(3,3,1),plot(x,spikemat_n,'color',[0.6 0.6 0.6])
                xlim([x(1) x(end)]);xlabel('Time(ms)');
                clear x
                hold on
                dot_raster_psth_plot(stimulus,lps,nreps,1,spktime_n,3,2,2,offset,SPK_Dur,randB,isonline);
                
                find_figure('spikes');
                x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
                subplot(3,3,4),plot(x,spikemat_n,'color',[0.6 0.6 0.6]);
                hold on
                cc=hsv(12);i=15;
                [~,c]=size(centers_n);
                for j=1:c
                    i=i-3;
                    subplot(3,3,4),plot(x,centers_n(:,j),'color',cc(i,:),'linewidth',2);
                    
                    hold on
                end
                xlim([x(1) x(end)]);xlabel('Time(ms)');
                clear x
                
                
                dot_raster_psth_plot(stimulus,lps,nreps,k_mean_n,spcounts_n,3,5,0,offset,SPK_Dur,randB,isonline);
                hold all
                
                %%%% for positive spikes%%%%%%%
                
                
                [nc_p,k_mean_p,k_centres_p]=nc_kmean(spikemat_p,n_spikes,p_spikes,2,13,3,2,8);
                [spcounts_p,spid_p,spid_unit_p,spmat_p,centers_p]=class_sorting(spikemat_p,spktime_p.cl1,spkid_p,nc_p,k_mean_p,2,2,8,k_centres_p);
                v = input('Do you want to do manual clustering (y/n)?','s');
                
                if strcmpi(v,'y')
                    spcounts_p=[]; centers_p=[]; spmat_p=[];
                    [spcounts_p,centers_p,spmat_p]=manual_cluster(spikemat_p,spkid_p,spktime_p.cl1);
                end
                
                
                find_figure('spikes');
                x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
                subplot(3,3,1),plot(x,spikemat_p,'color',[0.3 0.3 0.3])
                xlim([x(1) x(end)]);xlabel('Time(ms)');
                clear x
                hold on
                dot_raster_psth_plot(stimulus,lps,nreps,1,spktime_p,3,2,1,offset,SPK_Dur,randB,isonline);
                
                find_figure('spikes');
                x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
                subplot(3,3,7),plot(x,spikemat_p,'color',[0.3 0.3 0.3]);
                hold on
                
                cc=hsv(12);i=15;
                [~,c]=size(centers_p);
                for j=1:c
                    i=i-3;
                    subplot(3,3,7),plot(x,centers_p(:,j),'color',cc(i,:),'linewidth',2);
                    hold on
                end
                xlim([x(1) x(end)]);xlabel('Time(ms)');
                clear x
                hold off
                
                dot_raster_psth_plot(stimulus,lps,nreps,k_mean_p,spcounts_p,3,8,0,offset,SPK_Dur,randB,isonline);
                
                %%%%%%%%%%%end of positive+negative%%%%%%%%%%%%
            elseif strcmpi(w,'2')
                %%%%% only negative spikes%%%%%%%
                
                neg_spikes=1;
                pos_spikes=0;
                
                remove_spikes_n = input('Do you want to do remove  negative crossing spikes before clustering (y/n)?','s');
                
                if strcmpi(remove_spikes_n,'y')
                    idx=zeros(1,size(spikemat_n,2)); spikemat=[]; spkid=[]; spktime=[]; n_spikes =[];
                    spikemat = spikemat_n;
                    spkid = spkid_n;
                    spktime = spktime_n.cl1;
                    spikemat_n=[]; spkid_n=[]; spktime_n.cl1=[];
                    
                   [spikemat_n,spkid_n,spktime_n.cl1,~]= spike_removal_gui(spikemat,spkid,spktime,idx);
                   n_spikes = size(spikemat_n,2);
                end
                
                
                
                find_figure('Selection of nc');clf;
                find_figure('Selection of K for k-mean clustering');clf;
                [nc_n,k_mean_n,k_centres_n]=nc_kmean(spikemat_n,n_spikes,p_spikes,1,1,0,0,0);
                [spcounts_n,spid_n,spid_unit_n,spmat_n,centers_n]=class_sorting(spikemat_n,spktime_n.cl1,spkid_n,nc_n,k_mean_n,1,0,0,k_centres_n);
               
                
                find_figure('spikes');
                x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
                subplot(2,3,4),plot(x,spikemat_n,'color',[0.8 0.8 0.8]);
                hold on
                cc=hsv(12);i=12;
                [~,c]=size(centers_n);
                for j=1:c
                    subplot(2,3,4),plot(x,centers_n(:,j),'color',cc(i,:),'linewidth',2);
                    i=i-3;
                    hold on
                end
                xlim([x(1) x(end)]);xlabel('Time(ms)');
                clear x
                hold off
                
                dot_raster_psth_plot(stimulus,lps,nreps,k_mean_n,spcounts_n,2,5,0,offset,SPK_Dur,randB,isonline);
                
                %%%for manual clustering
                v = input('Do you want to do manual clustering (y/n)?','s');
                if strcmpi(v,'y')
                    spcounts_n=[]; centers_n=[]; spmat_n=[]; k_mean_n=[];
                    [spcounts_n,centers_n,spmat_n]= manual_cluster(spikemat_n,spkid_n,spktime_n.cl1);
                end
                find_figure('spikes');
                x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
                subplot(2,3,4),plot(x,spikemat_n,'color',[0.8 0.8 0.8]);
                hold on
                cc=hsv(12);i=12;
                [~,c]=size(centers_n);
                for j=1:c
                    subplot(2,3,4),plot(x,centers_n(:,j),'color',cc(i,:),'linewidth',2);
                    i=i-3;
                    hold on
                end
                xlim([x(1) x(end)]);xlabel('Time(ms)');
                clear x
                hold off
                h1=subplot(2,3,5);
                h2=subplot(2,3,6);
                cla(h1);cla(h2);
                k_mean_n=c;
                dot_raster_psth_plot(stimulus,lps,nreps,k_mean_n,spcounts_n,2,5,0,offset,SPK_Dur,randB,isonline);
                
                  %%% for manual clusters
            end
        else
            fprintf('\nNo spike crossing positive threshold ... \n')
            fprintf('\nDoing analysis for spike crossing negative threshold ... \n')
            
            remove_spikes_n = input('Do you want to do remove  negative crossing spikes before clustering (y/n)?','s');
                
                if strcmpi(remove_spikes_n,'y')
                    idx=zeros(1,size(spikemat_n,2)); spikemat=[]; spkid=[]; spktime=[]; n_spikes =[];
                    spikemat = spikemat_n;
                    spkid = spkid_n;
                    spktime = spktime_n.cl1;
                    spikemat_n=[]; spkid_n=[]; spktime_n.cl1=[];
                    
                   [spikemat_n,spkid_n,spktime_n.cl1,~]= spike_removal_gui(spikemat,spkid,spktime,idx);
                   n_spikes = size(spikemat_n,2);
                end
                
                
            find_figure('Selection of nc');clf;
            find_figure('Selection of K for k-mean clustering');clf;
            neg_spikes=1;
            pos_spikes=0;
              
            [nc_n,k_mean_n,k_centres_n]=nc_kmean(spikemat_n,n_spikes,p_spikes,1,1,0,0,0);
            [spcounts_n,spid_n,spid_unit_n,spmat_n,centers_n]=class_sorting(spikemat_n,spktime_n.cl1,spkid_n,nc_n,k_mean_n,1,0,0,k_centres_n);
            
            find_figure('spikes');
                x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
                subplot(2,3,4),plot(x,spikemat_n,'color',[0.8 0.8 0.8]);
                hold on
                cc=hsv(12);i=12;
                [~,c]=size(centers_n);
                for j=1:c
                    subplot(2,3,4),plot(x,centers_n(:,j),'color',cc(i,:),'linewidth',2);
                    i=i-3;
                    hold on
                end
                xlim([x(1) x(end)]);xlabel('Time(ms)');
                clear x
                hold off
                
                dot_raster_psth_plot(stimulus,lps,nreps,k_mean_n,spcounts_n,2,5,0,offset,SPK_Dur,randB,isonline);
                
            v = input('Do you want to do manual clustering (y/n)?','s');
            
            if strcmpi(v,'y')
                spcounts_n=[]; centers_n=[]; spmat_n=[]; k_mean_n=[];
                [spcounts_n,centers_n,spmat_n]=manual_cluster(spikemat_n,spkid_n,spktime_n.cl1);
            end
            find_figure('spikes');
            x=linspace(0,spike_len/rfs*1000,spike_len)-(stemplen*1000);
            subplot(2,3,4),plot(x,spikemat_n,'color',[0.8 0.8 0.8]);
            hold on
            cc=hsv(12);i=12;
            [~,c]=size(centers_n);
            for j=1:c
                subplot(2,3,4),plot(x,centers_n(:,j),'color',cc(i,:),'linewidth',2);
                i=i-3;
                hold on
            end
            xlim([x(1) x(end)]);xlabel('Time(ms)');
            clear x
            hold off
            k_mean_n=c;
            dot_raster_psth_plot(stimulus,lps,nreps,k_mean_n,spcounts_n,2,5,0,offset,SPK_Dur,randB,isonline);
        end
    end
end
% a=exist('w')
%               if a==0;
%                   neg_spikes=1;
%               end
