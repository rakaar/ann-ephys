function [nc,k_mean,k_centers_start]=nc_kmean(spikemat,n_spikes,p_spikes,row,dim1,dim2,dim3,dim4)


[nproj,total_spike]=size(spikemat);

spcorr=spikemat(:,:)*spikemat(:,:)';
[u,s,~]=svd(spcorr);

%%% calculating the no of principle components to be used%%%%%%%

find_figure('Selection of nc');

t=dim1;
for kk=1:4
    if kk>2
        t=kk+dim1+3;
    end
    subplot(row*2,6,t),plot(u(:,kk));
    t=t+1;
    title(sprintf('spike shape dim%i',kk))
end
clear t kk


%%% creating random spike noise


for jj=1:1000
    for kk=1:total_spike
        inx=randperm(nproj);
        rspkmat(:,kk)=spikemat(inx,kk);
    end
    rspcorr=rspkmat(:,:)*rspkmat(:,:)';
    [~,rs,~]=svd(rspcorr);
    rS(:,jj)=diag(rs);
end
clear jj kk
%%%
%%% plot eigen values

subplot(row,3,2+dim2),plot(mean(rS'),'r','linewidth',3);
hold on
subplot(row,3,2+dim2),plot(diag(s),'-o');
errorbar(mean(rS'),5*std(rS'));
title('plot of eigen values along with noise')

%%%%

p=diag(s);
sum_total=sum(p);
sum_n=0;
for i=1:nproj
    sum_n=sum(p(1:i));
    h(i)=sum_n/sum_total;
    subplot(row,3,3+dim2),plot(h,'-o');
    hold on
end
clear i sum_n

prompt='Select nc=';
nc=input(prompt);


%%%%Selection of nproj complete%%%%%%%%%%%

%%%%%%%%%%%%% figures for selection of K %%%%%%%%%%%%%%%%%%%%

if total_spike<10
    loops_n=total_spike;
else loops_n=10;
end

x=[]; y=[]; k_centers_start=[];

if loops_n==1 %%%clustering not possible
    k_mean=1;
    fprintf('Single spike present so no need for clustering \n')
    
    
else %% if spikes greater than 1
    
    proj_u=spikemat(:,:)'*u(:,1:nc);%% projection matrix using pca
    
    find_figure('Selection of K for k-mean clustering');
    
    if nc==1
        subplot(row,2,1+dim3),plot(proj_u(:,1),'.')
        title('cluster')
    elseif nc==2
        h2=subplot(row,2,1+dim3);
        plot(proj_u(:,1),proj_u(:,2),'.')
        title('cluster')
    elseif nc==3
        subplot(row,2,1+dim3),plot3(proj_u(:,1),proj_u(:,2),proj_u(:,3),'.')
        title('cluster')
    elseif nc==4
        dim=dim4;
        for ii=1:nc-2 %% to use nC3 combinations
            
            for jj=ii+1:nc-1
                if dim>1+dim4
                    dim=ii+jj+kk-4+dim4;
                end
                for kk=jj+1:nc
                    dim=dim+1;
                    eval(sprintf('t=proj_u(:,ii);'))
                    eval(sprintf('s=proj_u(:,jj);'))
                    eval(sprintf('r=proj_u(:,kk);'))
                    subplot(row*2,4,dim),plot3(t,s,r,'.')
                    %title(sprintf('nproj %i',dim))
                    
                end %% kk loop
            end %% jj loop
            
        end %% ii loop
        
    end
    
    
    
    for k_mean=2:loops_n
        clear idx c proj_class* pair_dist*
        [idx,c] = kmeans(proj_u, k_mean,'EmptyAction','singleton');%% k-mean clustering  'start','cluster',
        for kk=1:k_mean
            eval(sprintf('proj_class%i=proj_u(find(idx==%i),:);',kk,kk))%%%dividing the projection matrix
        end
        pair_sum=0;mean_pair=0;mean_c=0;
        for kk=1:k_mean
            eval(sprintf('pair_dist%i=mean(pdist(proj_class%i));',kk,kk))
            eval(sprintf('pair_dist=pair_dist%i;',kk))
            pair_dist(isnan(pair_dist))=0;
            pair_sum=pair_sum+pair_dist;
        end
        mean_pair=pair_sum/k_mean;
        mean_c=mean(pdist(c));
        inclass_incluster_ratio(k_mean)=mean_pair/mean_c;
    end
    inclass_incluster_ratio(isnan(inclass_incluster_ratio))=0;
    subplot(row,2,2+dim3),plot(2:loops_n,inclass_incluster_ratio(2:end))
    title('distance ratio for various cluster values')
    
    %%%%%
    
    
    if nc==4
        prompt='Select k for k-mean clustering=';
        k_mean=input(prompt);
        
    elseif nc==1
        prompt='Select k for k-mean clustering=';
        k_mean=input(prompt);
        
    elseif nc==3
        fprintf('select cluster centers along xy plane \n')
        figure
        g1=subplot(1,2,1);
        plot(proj_u(:,1),proj_u(:,2),'.')
        [x1,y]=getpts(g1);
        g2=subplot(1,2,2);
        plot(proj_u(:,1),proj_u(:,3),'.')
        hold on
        line([x1 x1],[min(proj_u(:,3)) max(proj_u(:,3))]);
        fprintf('select cluster centers along xz plane \n')
        [x2,z]=getpts(g2);
        
        ll=1;
        for j=1:length(x2)
            a = find(x1 >= x2(j)-0.005 & x1 <= x2(j)+0.005);
            if length(a)>1
                for l=1:length(a)
                    k_centers_start(ll,:)=[x1(a(l)) y(a(l)) z(j)];
                    ll=ll+1;
                end
            end
            k_centers_start(ll,:)=[x1(a) y(a) z(j)];
            ll=ll+1;
        end
        k_mean=size(k_centers_start,1);
        
    elseif nc==2
        fprintf('select cluster centers \n')
        find_figure('Selection of K for k-mean clustering');
        [x,y]=getpts(h2);
        k_centers_start=[x y];
        k_mean=size(k_centers_start,1);
    end
    
end
end




