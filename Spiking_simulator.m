%Script to implement the spiking network model found in
%%%
%Correlation-based model of artificially induced plasticity in motor cortex by a bilateral
%Brain-Machine Interface, Guillaume Lajoie, Nedialko Krouchev, John F. Kalaska, Adrienne Fairhall,
%Eberhard E. Fetz, under review, (2016), preprint: https://arxiv.org/abs/1609.00762
%%%

%Spike Time Dependent Plasticity (STDP) is implemented as an exponentially decaying rule based on pre-post
%spike-time difference and rate-based terms as well.

%This script simulates networks and compares to analytics

%--------------

%PARAMETERS

%simulation parameters
p.filename='spiking_test';
p.file_output='./';
p.function_path='./functions/'; %folder containing functions necessary for simulation
p.stim_switch=1; %Decides if conditioning is on (1) or off (0)
p.load_groups_switch=0; %THIS now controls an architecture load
p.load_path='./data.mat'; %must contain a structure "arch' with initial connectivity matrices J0,K0, J and K
p.T_total=4;%seconds ' length of total sim
p.T_batch_size=2;%seconds ' record things at every batch end
p.dt=0.5; %(ms) temporal rez
p.num_batch=p.T_total/p.T_batch_size;
p.num_batch_steps=p.T_batch_size*1000/p.dt;
p.ref_t=0; %refractory period (ms) (for spike generation.
addpath(p.function_path) %adding folder containing functions

%Stimulation switch
p.stim_del=0; %delay for stimulation (ms)
p.num_del_steps=p.stim_del/p.dt;
p.stim_size=10000; %Amount to rize rate of stimulated cells 
p.rec_neuron=1; %neuron to trigger stimulation
p.cols={'k','r','b'}; %colors for group-related plots

%number below needs to be divisible by 3 for groups 
p.N=60; %total number of neurons in net MUST BE DIVISIBLE BY 3

%Correlogram parameters (to ba sampled from spiking)
p.Plast_bin_center=-250:250; %resolution for building histograms of delta t's

%EXTERNAL RATE PARAMS
p.net_rate=3;%5; %(Hz) background rates for groups of neurons in network

%EXTERNAL RATE FUNCTIONS (SINE)
p.stim_type=0; %1=BOUTS 0=SINE
p.baseline=0; %(Hz)
p.ext_freq=[2,2,2]; %envelope frequencies (Hz);
p.ext_amplitude=[100,100,100]; %amplitude of envelope (Hz);
p.ext_offset=[0,2*pi/3,4*pi/3];

%EXTERNAL RATE FUNCTIONS (BOUTS)
% p.stim_type=1; %1=BOUTS 0=SINE
% p.baseline1=0; %(Hz) baseline firing rate --> Now see p.net_rate
% p.peaks1=[0.05,100]; %(Hz) minimum and maximum height of peaks
% p.peak_dur1=[5,100]; %(ms) minimum and maximum of peak width
% p.bout_dur1=[900,4000]; %(ms) minimum and maximum duration of activation bouts (ie duration of consecutive peaks)
% p.bout_freq1=0.13; %bouts / second (poisson distributed)
% p.ramp_dur=400; %ms
% p.ramp_ends=[0.1,5];

%SYNAPSES
p.tau_syn=5; %(ms) decay-time of synapse 
p.net_del=1; % (ms) mean of network synaptic delay
p.net_del_sig=1; %(ms) spread of network delays
p.dendrite_del=2; %(ms)
p.dendrite_ff_del=0; %(ms)

%STDP PARAMS
p.eta=5e-6;%5e-7; %learning rate time scale
p.w_in=0;%4;
p.w_out=0;%-0.5;
p.w_min=0;
p.w_max=0.02;%0.1;%0.02; %0.01
p.gamma=0.1;
% p.W=@STDPexp_del; %STDP rule
p.W=@(x,y)STDPJ(x,y,p.w_max,p.gamma);


%CONNECTIVITY MATRICES
p.net_connect_prob=.3; %prob of connection in network
p.net_w_init=0.005;%0.01;%0.015;%0.01, 0.015;

%GENERATING VARIABLES-------------------------------
%generating connectivity matrices
p.J0=double(rand(p.N)<p.net_connect_prob);
for i=1:p.N; p.J0(i,i)=0; end; %removing autapses
p.n_net=sum(sum(p.J0));
J=p.net_w_init*p.J0;

%_-_-_-_-_-_-_-_-_
% %Artificial assemblies
if p.load_groups_switch==1
    load(p.load_path);
    p.J0=arch.J0;
    J=arch.J;
end
%_-_-_-_-_-_-_-_-_

%generating delays
a=p.net_del-p.net_del_sig;
b=p.net_del+p.net_del_sig;
p.net_d=a + (b-a).*rand(p.N);
p.net_d=p.net_d.*p.J0;

%useful reference vector
ref_net=1:p.N;
p.group_ind={1:p.N/3,p.N/3+1:2*p.N/3,2*p.N/3+1:p.N};
p.back_buffer=p.Plast_bin_center(end); %time buffer for which to keep spikes from previous batches

%CONTAINERS
%Initialising containers -----------*
Jtraj=zeros(p.N,p.N,p.num_batch+1);
Jtraj(:,:,1)=J; %storing initial matrix
time=[0:p.num_batch]*p.T_batch_size; %time container
Rates=zeros(p.N,p.num_batch+1);
ISI_counts=zeros(3,3,length(p.Plast_bin_center),p.num_batch);
mean_ext_rates=zeros(3,p.num_batch);
RATES_array=zeros(3,p.num_batch_steps,p.num_batch);


%Synaptic spike queues and rolling containers
net_queu_len=round((max(max(p.net_d))+p.dendrite_ff_del)/p.dt+5);
net_syn_queue=p.T_batch_size*1000*10* ones(p.N,p.N,net_queu_len);
net_syn_entry=zeros(p.N,p.N,1);
net_syn_exit=zeros(p.N,p.N,1);
s=zeros(p.N,1); %synapses
r=zeros(p.N,1); %rates

%% Generating rates all at once
%container for full-length rates
num_total_t_steps=p.T_total*1000/p.dt;
Xtime=p.dt:p.dt:p.T_total*1000;
Xrates=zeros(3,num_total_t_steps);

%FOR SINE
if p.stim_type==0    
    for t=1:length(Xtime)
        x=ones(3,1)*t*p.dt/1000;
        Xrates(:,t)=(sin(2*pi*p.ext_freq'.*x+p.ext_offset')-4/5)*5;
        Xrates(:,t)=Xrates(:,t).*p.ext_amplitude';
        Xrates(:,t)=max(p.baseline*ones(3,1),Xrates(:,t));
    end 
end

%for BOUTS
if p.stim_type==1 
    for g=1:3 %loop over groups
        bout_starts=[];
        dur=0;
        t_stamp=0;
        tic
        while t_stamp<=p.T_total*1000 

            %drawing start of next bout
            delay1=exprnd(1000/p.bout_freq1);
            if isempty(bout_starts)
                bout_starts=[delay1];
            else
                bout_starts=[bout_starts bout_starts(end)+delay1+dur];
            end


            %drawing duration of next bout
            dur= p.bout_dur1(1) + (p.bout_dur1(2)-p.bout_dur1(1)).*rand;
            Sindex=find(Xtime>=bout_starts(end),1); %starting index
            Findex=Sindex; %finishing index
            GSindex=Sindex;
            while Findex*p.dt<=bout_starts(end)+dur
                peak_dur= p.peak_dur1(1) + (p.peak_dur1(2)-p.peak_dur1(1)).*rand;
                peak_height= p.peaks1(1) + (p.peaks1(2)-p.peaks1(1)).*rand;
                Findex=find(Xtime>=Xtime(Sindex)+peak_dur,1);
                Xrates(g,Sindex:Findex)=peak_height;
                Sindex=Findex+1;
            end
            t_stamp=bout_starts(end)+dur;

            %Working in ramps
            %up ramp
            ramp_index=GSindex:GSindex+p.ramp_dur/p.dt;
            Xrates(g,ramp_index)=(ramp_index-GSindex)*(p.ramp_ends(end)-p.ramp_ends(1))/(p.ramp_dur/p.dt)+p.ramp_ends(1);
            %down ramp
            ramp_index=Findex-p.ramp_dur/p.dt:Findex;
            Xrates(g,ramp_index)=(ramp_index-Findex)*(p.ramp_ends(1)-p.ramp_ends(end))/(p.ramp_dur/p.dt)+p.ramp_ends(1);
        end
        toc
    end
    Xrates=Xrates(:,1:num_total_t_steps);
    display('Done generating rates')
end

% rapackaging rates
for b=1:p.num_batch
    RATES_array(:,:,b)=Xrates(:,(b-1)*p.num_batch_steps+1:b*p.num_batch_steps);
    mean_ext_rates(:,b)=mean(Xrates,2);
end

%clearing rate-generation things
clear('Xrates','Xtime')
%% SIMULATION

for batch=1:p.num_batch
    display(['Processing batch # ' num2str(batch) '/' num2str(p.num_batch)])
    tic

    %-_-_-_-____-----_____-_-_-_
    %Pulling external rates from RATES_array
    %Generating random rates for this batch
    RATES=squeeze(RATES_array(:,:,batch));
    %containers for  BCI stim
    next_stim=[];
    %--------------

         
    %INITIATING LOCAL CONTAINERS-------
    %Spike arrays for network
    if batch>1
        past_spikes=spikes-p.T_batch_size*1000; %shifting spikes
        past_spike_entry=spike_entry;
    end
    spikes=zeros(p.N,100);%array of spikes with dealys added 
    spike_entry=ones(1,p.N);
    if batch>1
        for n=1:p.N
            spik=past_spikes(n,1:past_spike_entry(n)-1);
            spik=spik(spik>-p.back_buffer);
            spikes(n,1:length(spik))=spik;
            spike_entry(n)=length(spik)+1;
        end
    end
    
     %Synaptic spike queues and rolling containers
     if batch==1
        net_syn_queue=p.T_batch_size*1000*10* ones(p.N,p.N,net_queu_len);
        net_syn_entry=zeros(p.N,p.N,1);
        net_syn_exit=zeros(p.N,p.N,1);
     else
         net_syn_queue=net_syn_queue-p.T_batch_size*1000;
     end
    %----------------------------------
    
    %TEMPORARY___________________________
    r_temp=zeros(p.N,p.num_batch_steps);
    s_temp=zeros(p.N,p.num_batch_steps);
    %______________________________________
    
    %--------------
    %NET SIMULATION
    disp=0;
    for t=1:p.num_batch_steps
        t_now=t*p.dt;
        
        %progress display
        if 100*(t_now/1000/p.T_batch_size)-disp>10;
            disp=disp+10;
            display(['Processing batch # ' num2str(batch) '/' num2str(p.num_batch) '-->' num2str(disp) '% complete'])
        end
        
        %looping over net neurons and updating their weights according to
        %arriving spikes
        for n=1:p.N %loop over post-cells in net
            post_g=floor(3*(n-1)/p.N)+1; %group assignment
            %bianry spiking vectors
            net_spk_vect=zeros(p.N,1);
            
            %NET-------------
            pre_ind=ref_net(logical(p.J0(n,:))); %finding presyn neurons
            latest_spk=net_syn_queue(sub2ind(... %extracting latest spikes of pre cells
                size(net_syn_queue),...
                n*ones(size(pre_ind)),...
                pre_ind,...
                mod(net_syn_exit(n,pre_ind),net_queu_len)+1));
            spk_now_ind=pre_ind(latest_spk<=t_now); %index of pre cells that are spiking now
            if ~isempty(spk_now_ind)
                net_syn_queue(sub2ind(... %extracting latest spikes of pre cells
                size(net_syn_queue),...
                n*ones(size(spk_now_ind)),...
                spk_now_ind,...
                mod(net_syn_exit(n,spk_now_ind),net_queu_len)+1))...
                =p.T_batch_size*1000*10; %resetting old spikes
                net_spk_vect(spk_now_ind)=1; %storing spiking cells
                net_syn_exit(n,spk_now_ind)=net_syn_exit(n,spk_now_ind)+1;
            end
            %----------------
            
            %UPDATING SYNAPTIC VARIABLES
            s(n)=s(n)+J(n,:)*net_spk_vect/p.tau_syn;
            
            %$$$$$$$$$$$$$$$$$$$$$$$
            %PLASTICITY RATE-BASED
            %$$$$$$$$$$$$$$$$$$$$$$$
            %UPDATING J MATRIX with incoming spike rule
            if p.w_in~=0
                J(n,spk_now_ind)=J(n,spk_now_ind)+p.eta*p.w_in;
            end
            %$$$$$$$$$$$$$$$$$$$$$$$
             %$$$$$$$$$$$$$$$$$$$$$$$
             %$$$$$$$$$$$$$$$$$$$$$$$
             %$$$$$$$$$$$$$$$$$$$$$$$
             %STDP NOW !!!!!!!
            if spike_entry(n)>1 %check if current post-synaptic cell spiked at all in the past (otherwise there is no spike pair)
                %Post-synaptic spike train (all of it, adding dendritic delays)
                post_spk=spikes(n,1:spike_entry(n)-1)+p.dendrite_del; %Post spike arriving at synapse
                %$$$$$$$$$$$$$$$$$
                %DEPRESSION
                %STDP rule: running over PRESYN spikes from network
                for pre=spk_now_ind
                    pre_g=floor(3*(pre-1)/p.N)+1; %group assignment for current presynaptic cell
                    delta_t=t_now+p.dt/2-post_spk;
                    delta_t=delta_t(delta_t>0);
                    J(n,pre)=max(J(n,pre)+p.eta*sum(p.W(delta_t,J(n,pre))),p.w_min);
                    h=hist(delta_t,p.Plast_bin_center)';
                        h(1)=0;
                        h(end)=0;
                    ISI_counts(post_g,pre_g,:,batch)=squeeze(ISI_counts(post_g,pre_g,:,batch))+h;
                end
                %$$$$$$$$$$$$$$$$$
                
                %$$$$$$$$$$$$$$$$$
                %POTENTIATION
                %Checking if cell n (post) is spiking NOW after dendritic
                %delay
                spike_flag=0;
                postime=-1000;
                tik=length(post_spk);
                while tik>0 && post_spk(tik)>t_now-p.dt
                    if post_spk(tik)<=t_now && post_spk(tik)>t_now-p.dt
                        spike_flag=1;
                        postime=post_spk(tik);
                        break
                    end
                    tik=tik-1;
                end
                %Actual potentiation
                if spike_flag==1
                    for pre=pre_ind
                        pre_g=floor(3*(pre-1)/p.N)+1; %group assignment
                        del=p.net_d(n,pre);
                        pre_spk=spikes(pre,1:spike_entry(pre)-1)+del;
                        delta_t=pre_spk-postime; %creating delta t vector
                        delta_t=delta_t(delta_t<=0);
                        J(n,pre)=max(J(n,pre)+p.eta*sum(p.W(delta_t, J(n,pre))),p.w_min);
                        h=hist(delta_t,p.Plast_bin_center)';
                        h(1)=0;
                        h(end)=0;
                        ISI_counts(post_g,pre_g,:,batch)=squeeze(ISI_counts(post_g,pre_g,:,batch))+h;
                    end
                end
                %$$$$$$$$$$$$$$$$$
            end
            %$$$$$$$$$$$$$$$$$$$$$$$
            %$$$$$$$$$$$$$$$$$$$$$$$
            %$$$$$$$$$$$$$$$$$$$$$$$
        end %loop over post net neurons n
        
        %UPDATING NETWORK POISSON RATES---------------
        s=s-1/p.tau_syn*s*p.dt ; %synaptic decay
        ext=zeros(size(r)); %external rates
        ext(p.group_ind{1})=RATES(1,t)/1000;
        ext(p.group_ind{2})=RATES(2,t)/1000;
        ext(p.group_ind{3})=RATES(3,t)/1000;
        r=p.net_rate/1000+s+ext; %UPDATE
        %---------------------------------------------
        
        %Triggered Stim---------------------------------
        if p.stim_switch==1 && ~isempty(next_stim) && t_now>=next_stim(1)
            r(p.N/3+1:2*p.N/3)=r(p.N/3+1:2*p.N/3)+p.stim_size;
            if length(next_stim)>1
                next_stim=next_stim(2:end);
            else
                next_stim=[];
            end
        end
         %---------------------------------------------
        
        %ROLLING THE DIE TO SEE WHO SPIKES IN NET (and adding then to
        %queues)
        ind=ref_net(rand(p.N,1)<=r*p.dt); %indices of spiking cells
        if ~isempty(ind)
            
            %saving spike times now
            spikes(sub2ind(size(spikes),ind,spike_entry(ind)))=t_now+p.dt/2; %storing spikes
            spike_entry(ind)=spike_entry(ind)+1; %updating entry indices
            %resizing 'spikes' if too small
            if max(spike_entry)>size(spikes,2)
                spikes=[spikes,zeros(p.N,100)];
            end

            for n=ind %running over spiking cells now and storing spikes with according delays
                post_g=floor(3*(n-1)/p.N)+1; %group assignment
                
                Ntak=ref_net(logical(p.J0(:,n))); %post-syn indices
                Pretak=ref_net(logical(p.J0(n,:)));
                
                %UPDATING SYNAPTIC QUEUES
                net_syn_queue(sub2ind(size(net_syn_queue),Ntak,n*ones(size(Ntak)),...
                    mod(net_syn_entry(Ntak,n)',net_queu_len)+1))=t_now+p.dt/2+p.dendrite_ff_del+p.net_d(Ntak,n)'; %New here ! (added ff dendrite delay)
                net_syn_entry(Ntak,n)=net_syn_entry(Ntak,n)+1;
            end  %loop over spiking cells index

            %$$$$$$$$$$$$$$$$$$$$
            %PLASTICITY
            %UPDATING J MATRIX WITH OUTGOING SPIKE RULE
            if p.w_out~=0
                J(ind,:)=J(ind,:)+p.eta*p.w_out;
                J=J.*p.J0;
            end
            %$$$$$$$$$$$$$$$$$$$$
            
            %Chec rec cell for triggered stim------------------
            if ismember(p.rec_neuron,ind) %checking if REC neuron fires now
                next_stim=[next_stim t_now+p.dt/2+p.stim_del];
            end
            %------------------------------
        end
        
        %PLASTICITY $$$$$$$$
        % Apply J bounds
        J=min(J,p.w_max*ones(size(J)));
        J=max(J,p.w_min*ones(size(J)));
        %$$$$$$$$$$$$$$$$$$$$$
        
        %TEMPORARY STORING_______________
        r_temp(:,t)=r;
        s_temp(:,t)=s;
        %_____________________________________

        
    end %time loop for net simulation (batch)
    
    %Saving things from this batch
    Jtraj(:,:,batch+1)=J;
    %$$$--------------$$$
    
    %trimming spikes array
    spikes=spikes(:,1:max(spike_entry)-1);
    Rates(:,batch+1)=(spike_entry-1)/p.T_batch_size;
    
    toc %timer
    %saving things
    save([p.file_output, p.filename], 'p', 'Jtraj','time','Rates','spikes','batch','spike_entry','ISI_counts','mean_ext_rates','RATES_array')
end %batch loop


