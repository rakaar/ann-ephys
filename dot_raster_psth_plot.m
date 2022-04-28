function dot_raster_psth_plot(stimulus,lps,nreps,k_mean,spiketime,row,fig,spike_type,offset,SPK_Dur,randB,isonline)

% INPUT:stimulus=PP_PARAMS.protocol to get stimulus data
%       lps=lines per set i.e the no of frequency
%       nreps=no of repetion of same frequency
%       k_mean=to decide for negative or positive or cluster numbers
%       spiketime=input data which is spike time
%       offset_allow=is used to calculate offset in psth calculation
%       row=no of rows in which to plot the data either 2 or 3
%       fig=figure number for subplot
%       spike_type=is used to give different colours depending on positive
%                  or negative or cluster
%       offset=calculate the offset to be used in psth for various
%       frequency i.e. lps

if isfield(stimulus.stim_protocol,'level_hi')
    attn_gap=(stimulus.stim_protocol.level_hi-stimulus.stim_protocol.level_lo)/(stimulus.stim_protocol.num_att_steps-1);
end
stim_start=stimulus.stim_protocol.stim_start;
total_time=stimulus.stim_protocol.total_frames_lines;
binsz=0.01;
find_figure('spikes');
%%% dot raster plot
cc=hsv(12);i=15; %%% will be used to distinguish various cluster into colours

switch spike_type   %%1 for positive 2 for negative 0 for cluster
    case 1
        cc(15,:)=[0.3 0.3 0.3];
    case 2
        cc(15,:)=[0.6 0.6 0.6];
    case 0
        i=i-3;
end
for ii=1:k_mean
    for jj=1:lps
        for kk=1:nreps
            eval(sprintf('t=spiketime.cl%i.iter%i(:);',ii,(kk-1)*lps+jj))
            for ll=1:length(t)
                
                subplot(row,3,fig),hold on
                plot(t(ll),kk+(jj-1)*nreps,'.','color',cc(i,:));
                if isonline ~= 1
                    set(gca,'xlim',[stim_start-0.2 total_time]);
                end
                ylim([0 lps*nreps]);
                
            end
        end
    end
    i=i-3;
end


for jj=1:lps-1
    clear x y
    if rem(jj,2)==1
        if isonline== 1
            x=[0;0;SPK_Dur;SPK_Dur];
        else
            x=[stim_start-0.2;stim_start-0.2;total_time;total_time];
        end
        y=[nreps*jj;nreps+nreps*jj;nreps+nreps*jj;nreps*jj];
%         subplot(row,3,fig),hold on
%         h= patch(x,y,'r');
        
        ylim([0 lps*nreps]);
%         set(h,'facealpha',0.0000000000000000001);
    end
    
end
xlabel('Time(s)');
hold all
if isonline ~= 1
    stimuli_ht=0:nreps:lps*(nreps-1)+lps-1;
    stimulus_plot(stimulus,row,fig,lps*nreps,stimuli_ht,randB,lps*nreps);
end




%%% peristimuli time histogram

cc=hsv(12);i=15;

if stim_start-0.2>0 & isonline ~= 1
    edge_start=stim_start-0.2;
else
    edge_start=0;
end

if isonline == 1
    edges=edge_start:binsz:SPK_Dur;
else
    edges=edge_start:binsz:total_time; %%division of x axis into various bins
end


stimulus_height=zeros(1,k_mean);%% to be used to plot height of stimulus


if strcmpi(stimulus.type,'rsspip')
    
    if spike_type==1   %%1 for positive 2 for negative 0 for cluster
        cc(15,:)=[0.3 0.3 0.3];
    elseif spike_type==2
        cc(15,:)=[0.6 0.6 0.6];
    elseif spike_type==0
        i=i-3;
    end
    for ii=1:k_mean
        
        psth=zeros(lps,length(edges));
        psth_rss=zeros(1,length(edges));
        for jj=1:lps
            for kk=1:nreps
                eval(sprintf('xx=spiketime.cl%i.iter%i;',ii,(kk-1)*lps+jj))
                
                if ~isempty(xx)
                    psth(jj,:)=psth(jj,:)+histc(xx,edges);
                end
            end
            psth(jj,:)=smooth(psth(jj,:)/nreps*(1/binsz));
            
        end
        psth_rss=mean(psth);
        subplot(row,3,fig+1),plot(edges+binsz/2,psth_rss,'color',cc(i,:),'linewidth',2);
        
        stimulus_height(ii)=max(psth_rss)+.1*max(psth_rss);
        clear psth psth_rss
        i=i-3;
    end
else
    if spike_type==1   %%1 for positive 2 for negative 0 for cluster
        cc(15,:)=[0.3 0.3 0.3];
    elseif spike_type==2
        cc(15,:)=[0.6 0.6 0.6];
    elseif spike_type==0
        i=i-3;
    end
    for ii=1:k_mean
        psth=zeros(lps,length(edges));
        offset_psth=zeros(length(edges),lps);
        for jj=1:lps
            for kk=1:nreps
                eval(sprintf('xx=spiketime.cl%i.iter%i;',ii,(kk-1)*lps+jj))
                if ~isempty(xx)
                    psth(jj,:)=psth(jj,:)+histc(xx,edges);
                end
            end
            psth(jj,:)=smooth(psth(jj,:)/nreps*(1/binsz));
        end
        
        if isonline ~= 1
            psth=psth';
            offset_psth =  repmat(offset,length(edges),1) + psth ;
        else
            offset_psth=[];
            offset_psth=zeros(length(edges),2);
            offset_psth(:,1)=mean(psth(1:2:end,:))';
            offset_rand=max(offset_psth(:,1)) + 0.3 * max(offset_psth(:,1));
            offset_psth(:,2)=mean(psth(2:2:end,:))'+offset_rand;
            
        end
        
        subplot(row,3,fig+1),plot(edges+binsz/2,offset_psth,'color',cc(i,:),'linewidth',2);
        hold on
        stimulus_height(ii)=max(max(offset_psth))+0.1*max(max(offset_psth));
        clear psth offset_psth
        i=i-3;
    end
end
%bar(edges,psth');
stim_height=max(stimulus_height);
stim_height(isnan(stim_height))=1;
if stim_height==0
    stim_height=1;
end
if isonline ~= 1
    set(gca,'xlim',[stim_start-0.2 total_time]);
else
    set(gca,'xlim',[0 SPK_Dur]);
end
set(gca,'ylim',[0 stim_height]);
xlabel('Time(s)');
ylabel('sp/s');
hold all
if isonline ==1 %%%% only in case of random position of stimuli B we facemark the deviant
    x=[0;0;SPK_Dur;SPK_Dur];
    y=[offset_rand;stim_height;stim_height;offset_rand];
    subplot(row,3,fig+1),hold on
    h= patch(x,y,'r');
    set(h,'facealpha',0.0000000000000000005);
end

% if ~strcmpi(stimulus.type,'rsspip')&&lps~=1& isonline ~= 1  %% dividing the various frequencies in their respective repetitions
%     if rem(lps,2)==0
%         for jj=1:lps-2
%             clear x y
%             if rem(jj,2)==1
%                 if randB == 1
%                     x=[0;0;SPK_Dur;SPK_Dur];
%                 else
%                     x=[stim_start-0.2;stim_start-0.2;total_time;total_time];
%                 end
%                 y=[offset(jj+1);offset(jj+2);offset(jj+2);offset(jj+1)];
%                 subplot(row,3,fig+1),hold on
%                 h= patch(x,y,'r');
%                 set(h,'facealpha',0.05);
%             end
%         end
%         if randB == 1
%             x=[0;0;SPK_Dur;SPK_Dur];
%         else
%             x=[stim_start-0.2;stim_start-0.2;total_time;total_time];
%         end
%         y=[offset(end);stim_height;stim_height;offset(end)];
%         subplot(row,3,fig+1),hold on
%         h= patch(x,y,'r');
%         ylim([0 stim_height]);
%         set(h,'facealpha',0.05);
%     else
%         for jj=1:lps-1
%             clear x y
%             if rem(jj,2)==1
%                 if randB == 1
%                     x=[0;0;SPK_Dur;SPK_Dur];
%                 else
%                     x=[stim_start-0.2;stim_start-0.2;total_time;total_time];
%                 end
%                 y=[offset(jj+1);offset(jj+2);offset(jj+2);offset(jj+1)];
%                 hold on
%                 subplot(row,3,fig+1),hold on
%                 h= patch(x,y,'r');
%                 ylim([0 stim_height]);
%                 set(h,'facealpha',0.05);
%             end
%         end
%     end
% end
if ~strcmpi(stimulus.type,'rsspip')&& lps~=1 & isonline ~= 1  %% dividing the various frequencies in their respective repetitions
    if rem(lps,2)==0
        for jj=1:lps-2
            clear x y
            if rem(jj,2)==1
                x=[stim_start-0.2;stim_start-0.2;total_time;total_time];
                y=[offset(jj+1);offset(jj+2);offset(jj+2);offset(jj+1)];
                subplot(row,3,fig+1),hold on
                h= patch(x,y,'r');
                set(h,'facealpha',0.0000000000000000005);
            end
        end
        x=[stim_start-0.2;stim_start-0.2;total_time;total_time];
        y=[offset(end);stim_height;stim_height;offset(end)];
        subplot(row,3,fig+1),hold on
        h= patch(x,y,'r');
        ylim([0 stim_height]);
        set(h,'facealpha',0.0000000000000000005);
    else
        for jj=1:lps-1
            clear x y
            if rem(jj,2)==1
                x=[stim_start-0.2;stim_start-0.2;total_time;total_time];
                y=[offset(jj+1);offset(jj+2);offset(jj+2);offset(jj+1)];
                hold on
                subplot(row,3,fig+1),hold on
                h= patch(x,y,'r');
                ylim([0 stim_height]);
                set(h,'facealpha',0.0000000000000000005);
            end
        end
    end
end
hold all
if isonline ~= 1
    stimulus_plot(stimulus,row,fig+1,stim_height,offset,randB,lps*nreps);
end
end
