% script to find an equilibrium J given a C_hat and other parameters

%parameters
p.function_path='./functions/'; %folder containing functions necessary for simulation
addpath(p.function_path) %adding folder containing functions
p.N=60;
p.T_batch_size=2; %sec
p.stim_del=10; %delay for stimulation (ms)
p.net_del=2; %2 (ms) mean of network synaptic delay
p.net_del_sig=0; %(ms) spread of network delays
p.dendrite_del=2; %5 (ms)
p.dendrite_ff_del=2; %(ms)
p.w_min=0;
p.w_max=0.02;%0.1;%0.02; %0.01
p.net_connect_prob=.3; %prob of connection in network
p.net_w_init=(p.w_max-p.w_min)/2;%0.01;%0.015;%0.01, 0.015;
p.gamma=0.1; %learning rule param
p.eta=1e-6; %learning rate
p.expansion_order=1; %n
p.dt=0.5; %ms

%STDP parameters
%Gilson V
d.c_p=30;
d.tau_p=8.5;
d.c_d=20;
d.tau_d=17;

%For exp
% de.c_p=150;
% de.tau_p=20;
% de.c_d=200;
% de.tau_d=40;

%defining STDP rule
W_rule=@(delta_t,j)STDPJ(delta_t,j,p.w_max,p.gamma,d);
% W_rule=@(delta_t,j)STDPJ_exp(delta_t,j,p.w_max,p.gamma,de);


%% 
%Generating C_hat
dt=p.dt; %ms 
lags_xc=-300:p.dt:300;

%!!!!!!!!!!!!!!!!!!!!
%Gaussian-shaped C_hat
tau=50;
baseline_ext_xc=0.17;
peak_ext_xc=2.3;
mean_rate=sqrt(baseline_ext_xc/1000*p.dt)/p.dt*1000;
expected_count=mean_rate*p.T_batch_size;
% %---------
C_hat=zeros(3,3,length(lags_xc));
for pre=1:3
    for post=1:3
        if pre==post
            y=gaussmf(lags_xc,[tau,0])*peak_ext_xc/(1+baseline_ext_xc)+baseline_ext_xc;
%             y=baseline_ext_xc;
        else
            y=baseline_ext_xc;
        end
        C_hat(post,pre,:)=y;
    end
end 
%!!!!!!!!!!!!!!!!!!!!


%------


%%

%defining expansion rule
expansion=@(J,A)...
    C_hat2C(expected_count,C_hat,lags_xc,J,A,p.net_del+p.dendrite_ff_del,...
    p.net_del_sig,p.stim_del,p.expansion_order);

%Initial mean matrix
J_init=(p.w_max-p.w_min)/2*ones(3);
J=J_init*p.net_connect_prob*p.N/3;

%%
%NORMAL
%BBCI matrix
A=[0,0,0;0,0,0;0,0,0];
%Calling main function
display 'Computing Normal Equilibria'
[Jeq,C]=C_hat2Jeq(@(J)expansion(J,A),W_rule,J_init,lags_xc,p);

%BBCI

%BBCI matrix
A=[0,0,0;1,0,0;0,0,0];
%Calling main function
display 'Computing BBCI Equilibria'
[Jeq_dag,C_dag]=C_hat2Jeq(@(J)expansion(J,A),W_rule,J_init,lags_xc,p);

%% plots

%Cross-correlations
ma=max(max(max(C_dag)));
figure('Position', [100, 100, 1049, 895]);
for post=1:3
    for pre=1:3
        subplot(3,3,sub2ind([3,3],post,pre))
        set(gca,'FontSize',13)
        plot(lags_xc,squeeze(C_hat(post,pre,:)),'k')
        hold on
        plot(lags_xc,squeeze(C(post,pre,:)),'b')
        plot(lags_xc,squeeze(C_dag(post,pre,:)),'r')
        plot([0,0],[0,ma],'k--')
        hold off
        title([num2str(pre) '-->' num2str(post)])
        xlabel 'ISI(ms)'
        ylabel 'mean count'
        axis([lags_xc(1) lags_xc(end) 0 10])
    end
end
subplot(3,3,1)
legend('drive','normal','BBCI')
%%
%--------
%*************
%MISC
p.cols={'k','r','b'}; %colors for group-related plots
%input 9 colors to plot connections
p.group_cols=zeros(3,3,3);
p.group_cols(1,1,:)=[0,0,0]; %black for group 1
p.group_cols(2,2,:)=[1,0,0]; %red for group 2
p.group_cols(3,3,:)=[0,0,1]; %blue for group 3
p.group_cols(2,1,:)=[0.5,0,0];
p.group_cols(1,2,:)=[0.5,0.5,0];
p.group_cols(3,1,:)=[0,0,0.5];
p.group_cols(1,3,:)=[0,0.5,0.5];
p.group_cols(3,2,:)=[1,0,0.5];
p.group_cols(2,3,:)=[0.5,0,1];
%*************
%Equlibrium matrix
figure;
 %Plot mean matrix
subplot(2,2,1)
set(gca,'FontSize',13)
h=bar3(Jeq);
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
                         'facecolor',p.group_cols(ii,jj,:));
    end
end
rotate3d
%--
axis([0.5 3.5 0.5 3.5 p.w_min p.w_max])
xlabel 'pre'
ylabel 'post'
 title 'Before'
view(-30,60)

subplot(2,2,2)
set(gca,'FontSize',13)
h=bar3(Jeq_dag);
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
                         'facecolor',p.group_cols(ii,jj,:));
    end
end
rotate3d
%--
axis([0.5 3.5 0.5 3.5 p.w_min p.w_max])
xlabel 'pre'
ylabel 'post'
 title 'After'
view(-30,60)

% subplot(2,2,3)
% set(gca,'FontSize',10)
% title 'before'
% tix=cell(9,1);
% alpha={'a','b','c'};
% hold all
% counter=1;
% for pre=1:3
%     for post=1:3
%         bar(counter,Jeq(post,pre),'FaceColor',p.group_cols(post,pre,:));
%         tix{counter+1}=[alpha{pre} '-->' alpha{post}];
%         counter=counter+1;
%     end
% end
% hold off
% axis([0 10 p.w_min p.w_max])
% set(gca,'XTickLabel',tix)

% subplot(2,2,4)
% title 'after'
% set(gca,'FontSize',10)
% tix=cell(9,1);
% alpha={'a','b','c'};
% hold all
% counter=1;
% for pre=1:3
%     for post=1:3
%         bar(counter,Jeq_dag(post,pre),'FaceColor',p.group_cols(post,pre,:));
%         tix{counter+1}=[alpha{pre} '-->' alpha{post}];
%         counter=counter+1;
%     end
% end
% hold off
% axis([0 10 p.w_min p.w_max])
% set(gca,'XTickLabel',tix)

%changes
GA=(Jeq_dag-Jeq)/p.w_max;
subplot(2,2,[3,4])
tix=cell(9,1);
alpha={'a','b','c'};
hold all
counter=1;
for pre=1:3
    for post=1:3
        bar(counter,GA(post,pre),'FaceColor',p.group_cols(post,pre,:));
        tix{counter+1}=[];
        counter=counter+1;
    end
end
hold off
axis([0 10 -0.2 1])
set(gca,'XTickLabel',tix,'FontSize',13)
title('Normalized changes')

%% plot equilibrium dynamic curves
j_range=p.w_min:0.0001:p.w_max;
F=zeros(3,3,length(j_range));
F_dag=zeros(3,3,length(j_range));
for pre=1:3
    for post=1:3
        y=squeeze(C(post,pre,:));
        y_dag=squeeze(C_dag(post,pre,:));
        for j=1:length(j_range)
            F(post,pre,j)=sum(y'.*W_rule(lags_xc-p.dendrite_del+p.net_del,j_range(j)))*dt;
            F_dag(post,pre,j)=sum(y_dag'.*W_rule(lags_xc-p.dendrite_del+p.net_del,j_range(j)))*dt;
        end
    end
end

%plot of dynamic curves and equ changes
figure;

subplot(2,2,1)
title 'before'
set(gca,'FontSize',10)
tix=cell(9,1);
alpha={'a','b','c'};
hold all
for pre=1:3
    for post=1:3
        plot(j_range,squeeze(F(post,pre,:)),'Color',p.group_cols(post,pre,:))
        plot(Jeq(post,pre),0,'.','Color',p.group_cols(post,pre,:))
    end
end
plot(j_range,0*j_range,'k--')
hold off
set(gca,'XTickLabel',tix)

subplot(2,2,2)
title 'after'
set(gca,'FontSize',10)
tix=cell(9,1);
alpha={'a','b','c'};
hold all
for pre=1:3
    for post=1:3
        plot(j_range,squeeze(F_dag(post,pre,:)),'Color',p.group_cols(post,pre,:))
        plot(Jeq_dag(post,pre),0,'.','Color',p.group_cols(post,pre,:))
    end
end
plot(j_range,0*j_range,'k--')
hold off
set(gca,'XTickLabel',tix)

%changes
GA=(Jeq_dag-Jeq)/p.w_max;
subplot(2,2,[3,4])
set(gca,'FontSize',10)
tix=cell(9,1);
alpha={'a','b','c'};
hold all
counter=1;
for pre=1:3
    for post=1:3
        bar(counter,GA(post,pre),'FaceColor',p.group_cols(post,pre,:));
        tix{counter+1}=[alpha{pre} '-->' alpha{post}];
        counter=counter+1;
    end
end
hold off
% axis([0 10 -0.5 0.5])
set(gca,'XTickLabel',tix)




