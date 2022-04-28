function[offset]=psth_offset(spktime_n,spktime_p,binsz,stim_start,total_time,lps,nreps,isonline,SPK_Dur)

if stim_start-0.2>0 & isonline==0
    edge_start=stim_start-0.2;
else
    edge_start=0;
end

if isonline==1
   edges=edge_start:binsz:SPK_Dur;
else
   edges=edge_start:binsz:total_time;
end
psth_n=zeros(lps,length(edges));
psth_p=zeros(lps,length(edges));
offset=zeros(lps,1);
for jj=1:lps
    for kk=1:nreps
        eval(sprintf('xx=spktime_n.iter%i;',(kk-1)*lps+jj))
        eval(sprintf('yy=spktime_p.iter%i;',(kk-1)*lps+jj))
        if ~isempty(xx)
            psth_n(jj,:)=psth_n(jj,:)+histc(xx,edges);
        end
        if ~isempty(yy)
            psth_p(jj,:)=psth_p(jj,:)+histc(yy,edges);
        end
    end
    psth_n(jj,:)=smooth(psth_n(jj,:)/nreps*(1/binsz));
    psth_p(jj,:)=smooth(psth_p(jj,:)/nreps*(1/binsz));
end

offset_value=max(max(max(psth_n,psth_p)))+0.1*max(max(max(psth_n,psth_p)));

offset=offset_value*(0:lps-1);
end

