% Script to swipe over all delays and cross corr (C_hat) width for a given%
%Implements Gaussian-shaped C_hat

%Uses the functions: C_hat2C.m and C_hat2Jeq.m which do the expansion and
%convergnece to Equilibrium

%============================
%NETWORK PARAMETERS
%============================
%SIMULATION PARAMETERS
swipe_name='./test_analytic_swipe';
p.function_path='./functions/'; %folder containing functions necessary for simulation
addpath(p.function_path) %adding folder containing functions
save_flag=1;
p.T_batch_size=2;%seconds ' record things at every batch end
p.dt=0.5; %(ms) temporal rez
p.expansion_order=4; %order of cross-correlation expansion
%NETWORK SIZE
p.N=60;%60; %total number of neurons in net MUST BE DIVISIBLE BY 3
p.net_connect_prob=0.1;%0.10.3;
%SYNAPSES
p.tau_syn=5; %(ms) decay-time of synapse 
%STDP PARAMS
p.eta=5e-6;%5e-7; %learning rate time scale
p.w_in=0;%4;
p.w_out=0;%-0.5;
p.w_min=0;
p.w_max=0.0012;%0.02; %0.01
p.gamma=0.1;
p.net_w_init=p.w_max/4;%0.005;
p.net_del_sig=0; %(ms) spread of network delays
p.dendrite_ff_del=0; %(ms)
%============================
lags_xc=-300:p.dt:300; %time vector to evaluate XC's
baseline_ext_xc=0.17;
peak_ext_xc=2.3;
mean_rate=sqrt(baseline_ext_xc/1000*p.dt)/p.dt*1000;
expected_count=mean_rate*p.T_batch_size;
%============================
%============================

%============================
%SWIPE PARAMETERS
%============================
tau_xc_vect=[10,50];%[10:5:80];%[7,8,9,10,11];%[15 30 50 70 80 100]; %(ms) t-cst of external correlations
stim_del_vector=[15,10];%[0 1 5 10 20 50 80];%[0 1 5 10:10:250]; %stimulation delay vector
axonal_del_vect=[3];%[2,3,4]; %(ms) axonal delay
dendritic_back_del_vect=[3]; %(ms) back-propagation denrittic delay
%============================
%============================

%============================
%defining STDP rule and expansion functions
%============================
%STDP
d.c_p=30;
d.tau_p=8.5;
d.c_d=20;
d.tau_d=17;
W_rule=@(delta_t,j)STDPJ(delta_t,j,p.w_max,p.gamma,d);
% W_rule=@(delta_t,j)STDPJ_exp(delta_t,j,p.w_max,p.gamma,de);

%INITIAL MATRIX GUESS
J_init=(p.w_max-p.w_min)/2*ones(3);
J=J_init*p.net_connect_prob*p.N/3;
%============================
%============================

%============================
%CONTAINERS
Jeq_array=zeros(3,3,length(tau_xc_vect),length(stim_del_vector),length(axonal_del_vect),length(dendritic_back_del_vect));
Jeq_BCI_array=zeros(3,3,length(tau_xc_vect),length(stim_del_vector),length(axonal_del_vect),length(dendritic_back_del_vect));
Jeq_DIFF_array=zeros(3,3,length(tau_xc_vect),length(stim_del_vector),length(axonal_del_vect),length(dendritic_back_del_vect));
%============================
%============================


%% SWIPE
counter=0;
for t=1:length(tau_xc_vect)
    
    tau=tau_xc_vect(t);
    %------
    %Generating external cross-correlations
    C_hat=zeros(3,3,length(lags_xc));
    for pre=1:3
        for post=1:3
            if pre==post
                y=gaussmf(lags_xc,[tau,0])*peak_ext_xc/(1+baseline_ext_xc)+baseline_ext_xc;
            else
                y=baseline_ext_xc;
            end
            C_hat(post,pre,:)=y;
        end
    end 
    %------ 
    for d0=1:length(stim_del_vector)
        d_dag=stim_del_vector(d0);
        for d1=1:length(axonal_del_vect)
            da=axonal_del_vect(d1);
            for d2=1:length(dendritic_back_del_vect)
                dd_b=dendritic_back_del_vect(d2);
                counter=counter+1;

                %display progress
                tic
                display(['Processing ' num2str(counter) '/' num2str(prod(size(Jeq_array))/9)])

                %modifying parameter structure
                p.stim_del=d_dag;
                p.net_del=da; %2 (ms) mean of network synaptic delay
                p.dendrite_del=dd_b; %5 (ms)

                %EXPANSION HANDLE
                expansion=@(J,A)...
                    C_hat2C(expected_count,C_hat,lags_xc,J,A,p.net_del+p.dendrite_ff_del,...
                    p.net_del_sig,p.stim_del,p.expansion_order);

                %NORMAL EQ FINDER
                A=[0,0,0;0,0,0;0,0,0];
                [Jeq,C]=C_hat2Jeq(@(J)expansion(J,A),W_rule,J_init,lags_xc,p);
                Jeq_array(:,:,t,d0,d1,d2)=Jeq;

                %BBCI EQ FINDER
                A=[0,0,0;1,0,0;0,0,0];
                %Calling main function
                [Jeq_dag,C_dag]=C_hat2Jeq(@(J)expansion(J,A),W_rule,J_init,lags_xc,p);
                Jeq_BCI_array(:,:,t,d0,d1,d2)=Jeq_dag;

                %STRORING DIFFERENCE
                Jeq_DIFF_array(:,:,t,d0,d1,d2)=Jeq_dag-Jeq;

                toc
            end %d2
        end %d1
    end %d0
end %t

if save_flag==1
    save(swipe_name,'tau_xc_vect','stim_del_vector','axonal_del_vect','dendritic_back_del_vect','Jeq_array','Jeq_BCI_array','Jeq_DIFF_array','p','d','W_rule','baseline_ext_xc','peak_ext_xc','lags_xc');
end

%%
%==========
%PLOT PREP
para_names={'xc-width','stim del','axon del','dendrite bp del'};
group_cols=zeros(3,3,3);
group_cols(1,1,:)=[0,0,0]; %black for group 1
group_cols(2,2,:)=[1,0,0]; %red for group 2
group_cols(3,3,:)=[0,0,1]; %blue for group 3
group_cols(2,1,:)=[0.5,0,0];
group_cols(1,2,:)=[0.5,0.5,0];
group_cols(3,1,:)=[0,0,0.5];
group_cols(1,3,:)=[0,0.5,0.5];
group_cols(3,2,:)=[1,0,0.5];
group_cols(2,3,:)=[0.5,0,1];
%===========
%base choice index of all params
base_para=[1,1,1,1]; %[tau, d_dag, d_axonal, d_den_back, d_den_ff]
base_indexing={':',':',base_para(1),base_para(2),base_para(3),base_para(4)};
centralized_para_vects={tau_xc_vect,stim_del_vector,axonal_del_vect,dendritic_back_del_vect};

% Plotting the difference between connections for one chosen parameter
%Choice of paramter to plot
para_choice=2; %1=tau_xc, 2=stim delay, 3=axonal delay, 4=dendritic back prop, 5=dendritic feed forward
pull_indexing=base_indexing;
pull_indexing{para_choice+2}=':';

%Pulling data
D=squeeze(Jeq_DIFF_array(pull_indexing{:}))/(p.w_max-p.w_min);

figure;
set(gca,'FontSize',13)
hold all
tix=cell(9,1);
alpha={'a','b','c'};
counter=1;
for pre=1:3
    for post=1:3
        plot(centralized_para_vects{para_choice}, squeeze(D(post,pre,:)),'Color',group_cols(post,pre,:),...
            'LineWidth',2)
        tix{counter}=[alpha{pre} '-->' alpha{post}];
        counter=counter+1;
    end
end
hold off
xlabel([para_names{para_choice} ' (ms)'])
ylabel 'Relative difference'
% legend(tix);
K=1:4;
K=K(K~=para_choice);
tit=[];
for k=K
    tit=[tit para_names{k} '=' num2str(centralized_para_vects{k}(base_para(k))) ' (ms), '];
end
title(tit)

%% Surface plot for chosen pre/post
% Surface Plot of the difference between connections for two chosen parameter, for a choice of synapse
%Choice of paramter to plot
pre=1;
post=2;
para_choice1=1; %1=tau_xc, 2=stim delay, 3=axonal delay, 4=dendritic back prop
para_choice2=2;
pull_indexing=base_indexing;
%fixing indexing
pull_indexing{para_choice1+2}=':';
pull_indexing{para_choice2+2}=':';
pull_indexing{1}=post;
pull_indexing{2}=pre;

%Pulling data
D=squeeze(Jeq_DIFF_array(pull_indexing{:}))/(p.w_max-p.w_min);
x=centralized_para_vects{para_choice2};
y=centralized_para_vects{para_choice1};
figure;
set(gca,'FontSize',13)
surf(x,y,D)
shading interp
xlabel([para_names{para_choice2} ' (ms)'])
ylabel([para_names{para_choice1} ' (ms)'])
K=1:4;
K=K(K~=para_choice1);
K=K(K~=para_choice2);
tit=[num2str(pre) '-->' num2str(post) ', '];
for k=K
    tit=[tit para_names{k} '=' num2str(centralized_para_vects{k}(base_para(k))) ' (ms), '];
end
title(tit)
axis([x(1) x(end) y(1) y(end) -0.5 1])
caxis([-0.2,1])
colorbar

%% Grid color plot
% Color Plot of the difference between connections for two chosen
% parameter, for all synaptic pairs

base_para=[1,1,1,1]; %[tau, d_dag, d_axonal, d_den_back, d_den_ff]
base_indexing={':',':',base_para(1),base_para(2),base_para(3),base_para(4)};

para_choice1=1; %1=tau_xc, 2=stim delay, 3=axonal delay, 4=dendritic back prop
para_choice2=2;

figure;
for post=1:3
    for pre=1:3;

        pull_indexing=base_indexing;
        %fixing indexing
        pull_indexing{para_choice1+2}=':';
        pull_indexing{para_choice2+2}=':';
        pull_indexing{1}=post;
        pull_indexing{2}=pre;

        %Pulling data
        D=squeeze(Jeq_DIFF_array(pull_indexing{:}))'/(p.w_max-p.w_min);

        subplot(3,3,sub2ind([3,3],pre,post))
        set(gca,'FontSize',13)
        x=centralized_para_vects{para_choice2};
        y=centralized_para_vects{para_choice1};
        surf(x,y,D')
        shading interp
%         xlabel([para_names{para_choice1} ' (ms)'])
%         ylabel([para_names{para_choice2} ' (ms)'])
        caxis([-0.2,0.6])
        axis([x(1) x(end) y(1) y(end) -0.5 1])
        K=1:4;
        K=K(K~=para_choice1);
        K=K(K~=para_choice2);
        tit=[num2str(pre) '-->' num2str(post)];
%         for k=K
%             tit=[tit para_names{k} '=' num2str(centralized_para_vects{k}(base_para(k))) ' (ms), '];
%         end
%         title(tit)
% axis off
    end
end
subplot(3,3,1)
tit=[num2str(1) '-->' num2str(1) ', '];
for k=K
    tit=[tit para_names{k} '=' num2str(centralized_para_vects{k}(base_para(k))) ' (ms), '];
end
% title(tit)

%%  Plot of bar graph for before and after for chosen parameter point
base_para=[1,1,1,1]; %[tau, d_dag, d_axonal, d_den_back, d_den_ff]
base_indexing={':',':',base_para(1),base_para(2),base_para(3),base_para(4)};

%pulling data
D=squeeze(Jeq_DIFF_array(base_indexing{:}))/(p.w_max-p.w_min);
Be=squeeze(Jeq_array(base_indexing{:}));
Af=squeeze(Jeq_BCI_array(base_indexing{:}));

%Potting cross-correlations
tau=centralized_para_vects{1}(base_para(1));
%------
%Generating external cross-correlations
C_hat=zeros(3,3,length(lags_xc));
for pre=1:3
    for post=1:3
        if pre==post
            y=gaussmf(lags_xc,[tau,0])*peak_ext_xc/(1+baseline_ext_xc)+baseline_ext_xc;
        else
            y=baseline_ext_xc;
        end
        C_hat(post,pre,:)=y;
    end
end 
%------ 
figure;
ma=max(max(max(C_hat)));
for post=1:3
    for pre=1:3
        subplot(3,3,sub2ind([3,3],post,pre))
        set(gca,'FontSize',13)
        plot(lags_xc,squeeze(C_hat(post,pre,:)),'k')
        hold on
        plot([0,0],[0,ma],'k--')
        hold off
%         title([num2str(pre) '-->' num2str(post)])
%         xlabel 'ISI(ms)'
%         ylabel 'mean count'
%         axis([lags_xc(1) lags_xc(end) 0 ma])
axis off
    end
end


%Before and After connectivity matrix
figure;
subplot(2,2,[1,2])
set(gca,'FontSize',15)
tix=cell(9,1);
alpha={'a','b','c'};
hold all
counter=1;
for pre=1:3
    for post=1:3
        bar(counter,D(post,pre),'FaceColor',squeeze(group_cols(post,pre,:)));
        tix{counter}=[alpha{pre} '-->' alpha{post}];
        counter=counter+1;
    end
end
hold off
set(gca,'XTickLabel',tix)
axis([1 9 -0.5 1])
K=1:4;
tit=[];
for k=K
    tit=[tit para_names{k} '=' num2str(centralized_para_vects{k}(base_para(k))) ' (ms), '];
end
title(tit)

subplot(2,2,3)
set(gca,'FontSize',15)
h=bar3(Be);
%-- Group Colors
cnt = 0;
for jj = 1:length(h)
    xd = get(h(jj),'xdata');
    yd = get(h(jj),'ydata');
    zd = get(h(jj),'zdata');
    delete(h(jj))    
    idx = [0;find(all(isnan(xd),2))];
    if jj == 1
        S = zeros(length(h)*(length(idx)-1),1);
    end
    for ii = 1:length(idx)-1
        cnt = cnt + 1;
        S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
                         yd(idx(ii)+1:idx(ii+1)-1,:),...
                         zd(idx(ii)+1:idx(ii+1)-1,:),...
                         'facecolor',group_cols(ii,jj,:));
    end
end
rotate3d
%--
axis([0.5 3.5 0.5 3.5 p.w_min p.w_max])
xlabel 'pre'
ylabel 'post'
title 'Before Stimulation'
view(-30,60)

subplot(2,2,4)
set(gca,'FontSize',15)
h=bar3(Af);
%-- Group Colors
cnt = 0;
for jj = 1:length(h)
    xd = get(h(jj),'xdata');
    yd = get(h(jj),'ydata');
    zd = get(h(jj),'zdata');
    delete(h(jj))    
    idx = [0;find(all(isnan(xd),2))];
    if jj == 1
        S = zeros(length(h)*(length(idx)-1),1);
    end
    for ii = 1:length(idx)-1
        cnt = cnt + 1;
        S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
                         yd(idx(ii)+1:idx(ii+1)-1,:),...
                         zd(idx(ii)+1:idx(ii+1)-1,:),...
                         'facecolor',group_cols(ii,jj,:));
    end
end
rotate3d
%--
axis([0.5 3.5 0.5 3.5 p.w_min p.w_max])
xlabel 'pre'
ylabel 'post'
title 'After Stimulation'
view(-30,60)


