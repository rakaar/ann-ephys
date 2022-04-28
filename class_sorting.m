function [spcounts,spid,spid_unit,spmat,centers]=class_sorting(spikemat,spktime,spkid,nc,k_mean,row,dim1,dim2,k_centers_start)
%%%% principle component analysis done on spikemat(consisting of 'nproj' main
%%%% component)
%%%%k-mean clustering(where k=nc) is done on pca matrix proj_u
%%%fprintf('\ Doing Spike Sorting in classes... \n')
fprintf('Doing Spike Sorting into class ... \n')
[~,total_spike]=size(spikemat);
iters=length(fieldnames(spktime));%% to calculate number of iterations
spcorr=spikemat(:,:)*spikemat(:,:)';
[u,~,~]=svd(spcorr);
proj_u=spikemat(:,:)'*u(:,1:nc);

if nc==2 | nc==3 %%% to start clustering from a given centers
    [idx,c]=kmeans(proj_u,k_mean,'EmptyAction','singleton','start',k_centers_start);
else
    [idx,c]=kmeans(proj_u,k_mean,'EmptyAction','singleton');
end

centers=u(:,1:nc)*c';%%%gives the mean spike of various class


%%%%%%dividing the spike matrix and projection matrix into the classes
%%%%%%defined by k-mean clustering

for kk=1:k_mean
    eval(sprintf('proj_u_class%i=proj_u(find(idx==%i),:);',kk,kk))%%%dividing the projection matrix
end
for kk=1:k_mean
    eval(sprintf('class_set%i=find(idx==%i);',kk,kk))
end




find_figure('Selection of K for k-mean clustering');

if k_mean~=1
    
    cc=hsv(12);i=12; %%% will be used to distinguish various cluster into colours
    switch nc
        
        case 1
            for kk=1:k_mean
                eval(sprintf('t=proj_u_class%i(:,1);',kk))
                subplot(row,2,1+dim1),plot(t,'.','color',cc(i,:));
                i=i-2;
                hold on
            end
            
            
        case 2
            
            for ii=1:nc-1 %% to use nC2 combinations
                for jj=ii+1:nc
                    i=12;
                    for kk=1:k_mean
                        eval(sprintf('t=proj_u_class%i(:,ii);',kk))
                        eval(sprintf('s=proj_u_class%i(:,jj);',kk))
                        
                        subplot(row,2,1+dim1),plot(t,s,'.','color',cc(i,:));
                        
                        i=i-3;
                        hold on
                    end %% nc loop
                    hold off
                end %% jj loop
            end %% ii loop
            
            
        case 3
            
            for ii=1:nc-1 %% to use nC2 combinations
                for jj=ii+1:nc
                    i=12;
                    for ll=jj+1:nc
                        for kk=1:k_mean
                            eval(sprintf('t=proj_u_class%i(:,ii);',kk))
                            eval(sprintf('s=proj_u_class%i(:,jj);',kk))
                            eval(sprintf('r=proj_u_class%i(:,ll);',kk))
                            subplot(row,2,1+dim1),plot3(t,s,r,'.','color',cc(i,:));
                            
                            i=i-3;
                            hold on
                        end %% nc loop
                    end %% ll loop
                end %% jj loop
            end %% ii loop
            
            
        case 4
            
            dim=dim2;
            for ii=1:nc-2 %% to use nC3 combinations
                for jj=ii+1:nc-1
                    if dim>1+dim2
                        dim=ii+jj+ll+dim2-4;
                    end
                    for ll=jj+1:nc
                        dim=dim+1;
                        i=12;
                        for kk=1:k_mean
                            eval(sprintf('t=proj_u_class%i(:,ii);',kk))
                            eval(sprintf('s=proj_u_class%i(:,jj);',kk))
                            eval(sprintf('r=proj_u_class%i(:,ll);',kk))
                            subplot(row*2,4,dim),plot3(t,s,r,'.','color',cc(i,:));
                            %%title(sprintf('nproj %i',nproj))
                            i=i-3;
                            hold on
                        end %% nc loop
                        
                    end %% ll loop
                end %% jj loop
            end %%ii loop
            
    end
end



%%%%%%%%%%%dividing spike time into various classes
%%%% cannot be done directly so using spkid containing spike no of each
%%%% iteration is used to divide into classes
for kk=1:iters
    eval(sprintf('itersids=spkid.iter%i(:);',kk))
    for sn=1:length(itersids)
        if k_mean==1
            find(class_set1==itersids(sn));
            iterclass(sn)=1;
            
        elseif k_mean==2
            if find(class_set1==itersids(sn))
                iterclass(sn)=1;
            elseif find(class_set2==itersids(sn))
                iterclass(sn)=2;
            end
            
        elseif k_mean==3
            if find(class_set1==itersids(sn))
                iterclass(sn)=1;
            elseif find(class_set2==itersids(sn))
                iterclass(sn)=2;
            elseif find(class_set3==itersids(sn))
                iterclass(sn)=3;
            end
            
        elseif k_mean==4
            if find(class_set1==itersids(sn))
                iterclass(sn)=1;
            elseif find(class_set2==itersids(sn))
                iterclass(sn)=2;
            elseif find(class_set3==itersids(sn))
                iterclass(sn)=3;
            elseif find(class_set4==itersids(sn))
                iterclass(sn)=4;
            end
            
        end
        
    end
    
    if exist('iterclass')
        eval(sprintf('spkiterclass.iter%i=iterclass;',kk))
        clear iterclass
    else
        eval(sprintf('spkiterclass.iter%i=[];',kk))
    end
end
%%%%% dividing the spike time instant into various classes using the
%%%%% spkiterclass
for unit=1:k_mean
    eval(sprintf('spid_unit%i=[];',unit))
    for kk=1:iters
        eval(sprintf('spcounts.cl%i.iter%i=spktime.iter%i(find(spkiterclass.iter%i==unit));',unit,kk,kk,kk));
        eval(sprintf('spid.cl%i.iter%i=spkid.iter%i(find(spkiterclass.iter%i==unit));',unit,kk,kk,kk))
        eval(sprintf('spid_unit%i=[spid_unit%i spid.cl%i.iter%i];',unit,unit,unit,kk))
        eval(sprintf('spid_unit.cl%i=spid_unit%i;',unit,unit))
    end
    eval(sprintf('[xx yy]=size(spid_unit.cl%i);',unit));
    for tt=1:yy
        eval(sprintf('col_spikemat=spid_unit.cl%i(%i);',unit,tt))
        eval(sprintf('spmat.cl%i(:,%i)=spikemat(:,%i);',unit,tt,col_spikemat))
    end
end
end


