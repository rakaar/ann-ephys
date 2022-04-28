function[]=stimulus_plot(stimulus,row,spnum,stim_height,offset,randB,iters)

%INPUT: stimulus=PP_PARAMS.protocol to get stimulus data
%       row=to decide the no of rows in subplot either 2 for only negative
%           or 3 when both negative and positive
%       spnum=to give the figure no where to plot the stimulus
%       stim_height=the height of the stimulus is provided for dot raster
%                   its lps*nreps and for psth its decided by max value of offset_psth
%In this function the time instant of the stimulus being played is ploted
%on the dot raster as well as psth plot
%

find_figure('spikes');

switch stimulus.type
    
    case 'NoisePip'
        stim_start=stimulus.stim_protocol.stim_start;
        npipdur=stimulus.stim_protocol.pip_dur;
        ipi=stimulus.stim_protocol.ipi_val;
        npips=stimulus.stim_protocol.npips_per_train;
        %stim_length=pip_dur+(npips_per_train-1)*ipi_val;
        
        for jj=1:npips
            pip_start=stim_start+(npipdur+ipi)*(jj-1);
            pip_stop=stim_start+(npipdur+ipi)*(jj-1)+npipdur;
            x=[pip_start;pip_start;pip_stop;pip_stop];
            y=[0;stim_height;stim_height;0];
            subplot(row,3,spnum),hold on
            h1=patch(x,y,'b');
            set(h1,'facealpha',0.01)
            
        end
        
    case 'TonePipSweep'
        stim_start=stimulus.stim_protocol.stim_start;
        tpipdur=stimulus.stim_protocol.pip_dur;
        ipi=stimulus.stim_protocol.ipi_val;
        npips=stimulus.stim_protocol.npips_per_train;
        %stim_length=pip_dur+(npips_per_train-1)*ipi_val;
        
        for jj=1:npips
            pip_start=stim_start+(tpipdur+ipi)*(jj-1);
            pip_stop=stim_start+(tpipdur+ipi)*(jj-1)+tpipdur;
            x=[pip_start;pip_start;pip_stop;pip_stop];
            y=[0;stim_height;stim_height;0];
            subplot(row,3,spnum),hold on
            h1=patch(x,y,'b');
            set(h1,'facealpha',0.01)
            
        end
        
    case 'SSA' %%% assumes A and B tokens of same length
        stim_start=stimulus.stim_protocol.stim_start;
        if isfield(stimulus.stim_protocol,'rand_nstim_check')
            var_stim=stimulus.stim_protocol.rand_nstim_check;
        else
            var_stim=0;
        end
        
        if randB~=1 & var_stim~=1
            
            for jj=1:stimulus.stim_protocol.total_stims
                token_start=stim_start+(stimulus.stim_protocol.gap_dur+stimulus.stim_protocol.stim_Adur)*(jj-1);
                token_stop=stim_start+(stimulus.stim_protocol.gap_dur+stimulus.stim_protocol.stim_Adur)*(jj-1)+stimulus.stim_protocol.stim_Adur;
                if jj==stimulus.stim_protocol.stimB_loc
                    x=[token_start;token_start;token_stop;token_stop];
                    y=[0;stim_height;stim_height;0];
                    subplot(row,3,spnum),hold on
                    h1=patch(x,y,'r');
                    set(h1,'facealpha',0.001)
                    
                else
                    x=[token_start;token_start;token_stop;token_stop];
                    y=[0;stim_height;stim_height;0];
                    subplot(row,3,spnum),hold on
                    h1=patch(x,y,'b');
                    set(h1,'facealpha',0.001)
                    
                end
            end
            
        else
            
            if randB==1
                stimB_loc = stimulus.stim_protocol.rand_stimB_loc;
            else
                stimB_loc = repmat(stimulus.stim_protocol.stimB_loc,1,iters);
            end
            
            if var_stim==1
                stimulus_nos=stimulus.stim_protocol.rand_nstim_num;
            else
                stimulus_nos= repmat(stimulus.stim_protocol.total_stims,1,iters);
            end
            
            
            for i=1:length(stimulus_nos)
                total_stims=0;
                eval(sprintf('total_stims=stimulus_nos(%i);',i))
                eval(sprintf('stimB_pos=stimB_loc(%i);',i))
                
                for j=1:total_stims
                    token_start=stim_start+(stimulus.stim_protocol.gap_dur+stimulus.stim_protocol.stim_Adur)*(j-1);
                    token_stop=stim_start+(stimulus.stim_protocol.gap_dur+stimulus.stim_protocol.stim_Adur)*(j-1)+stimulus.stim_protocol.stim_Adur;
                    
                    x=[token_start;token_start;token_stop;token_stop];
                    eval(sprintf('y_height_low=offset(%i);',i))
                    
                    if i ~= length(stimulus_nos)
                        eval(sprintf('y_height_high=offset(%i);',i+1))
                    else
                        y_height_high=stim_height;
                    end
                    
                    y=[y_height_low;y_height_high;y_height_high;y_height_low];
                    subplot(row,3,spnum),hold on
                    
                    if j== stimB_pos
                        h1=patch(x,y,'r');
                        set(h1,'facealpha',0.3)
                    else
                        h1=patch(x,y,'b');
                        set(h1,'facealpha',0.1)
                    end
                end
            end
            
        end
        
    case 'RSSPip'
        stim_start=stimulus.stim_protocol.stim_start;
        npipdur=stimulus.stim_protocol.pip_dur;
        ipi=stimulus.stim_protocol.ipi_val;
        npips=stimulus.stim_protocol.npips_per_train;
        %stim_length=pip_dur+(npips_per_train-1)*ipi_val;
        
        for jj=1:npips
            pip_start=stim_start+(npipdur+ipi)*(jj-1);
            pip_stop=stim_start+(npipdur+ipi)*(jj-1)+npipdur;
            x=[pip_start;pip_start;pip_stop;pip_stop];
            y=[0;stim_height;stim_height;0];
            subplot(row,3,spnum),hold on
            h1=patch(x,y,'b');
            set(h1,'facealpha',0.1)
            
        end
        
end
end