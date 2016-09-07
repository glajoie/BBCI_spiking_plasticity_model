%GL Feb 2016

%Function that does expansion of estimate for network cross-corr C(u) to
%desired order

%INPUTS:
%expected_count= expected number of spikes in epoch from each group.
%C_hat: external cross-corr 
%lags_xc: time axis for C_hat
%J: connectivity matrix average TOTAL INCOMING SYNAPSES
%A: spike triggered stimulation matrix
%da: mean axonal delay + ff dendritic delay (if present)
%da_sig: sigma of axonal delay 
%ds: delay of spike-triggered stimulation
%n: order of exnpansion

%OUTPUT
%C: network cross-correlation structure, defined with time-dicretization
%given by lags_xc

%============================

function [C]=C_hat2C(expected_count,C_hat,lags_xc,J,A,da,da_sig,ds,n)

%n is the order of expansion
dt=lags_xc(2)-lags_xc(1);

%% New code
C=C_hat;
normal_mats={J,A};
transpose_mats={J',A'};
for k=1:n
    k_combos=allVL1(4,k); %all possible pair of integers (i,j,x,y) summing to k
    for k_runner=1:size(k_combos,1) %running over all combos of i and j summin to k
        
        %extracting current values of i and j
        i=k_combos(k_runner,1);
        j=k_combos(k_runner,2);
        x=k_combos(k_runner,3);
        y=k_combos(k_runner,4);
        
        %extracting shifted C_hat
        t_shift=(i-j)*da+(x-y)*ds; %current time-shift
        C_temp=dshift(C_hat,t_shift,dt); %shifted copy of C_hat
        
       %creating left indices
        left_permut=[];
        if i+x>0
            temp=zeros(1,i+x);
            [temp(1:i)]=1; %J
            [temp(i+1:i+x)]=2; %A
            left_permut=unique(perms(temp),'rows');
        end
        
        %creating right indices
        right_permut=[];
        if j+y>0
            temp=zeros(1,j+y);
            [temp(1:j)]=1; %J'
            [temp(j+1:j+y)]=2; %A'
            right_permut=unique(perms(temp),'rows');
        end
        
        %Adding up terms
        for left=1:nchoosek(i+x,i) %loop over left permutations
            
            %extracting current left term
            Gleft=eye(size(J));
            if ~isempty(left_permut)
                gf=normal_mats(left_permut(left,:)); 
                for tak=1:length(gf) %multiply current permutation
                    Gleft=Gleft*gf{tak};
                end
            end

            for right=1:nchoosek(j+y,j) %loop over right permutations
                
                 %extracting current left term
                Gright=eye(size(J));
                if ~isempty(right_permut)
                    gf=transpose_mats(right_permut(right,:)); 
                    for tak=1:length(gf) %multiply current permutation
                        Gright=Gright*gf{tak};
                    end
                end
                
                %Adding up to expansion !
                for u=1:size(C,3);  %running over time index and making matrix product
                    C(:,:,u)=C(:,:,u)+Gleft*squeeze(C_temp(:,:,u))*Gright;
                end%u
                
            end%right
        end%left
    end %k_runner
end %k


%% Delta-function corrections self-consistent
if sum(sum(A~=0))>0 %do delta-fct adjustments
    
    %creating a delta-only corr matix
    C_delta=0*C_hat;
    [post,pre]=ind2sub(size(A),find(A~=0));
    %autocorr: %adding delta peak with width due to distribution of delays
    %in stimulated areas
    u0=(length(lags_xc)-1)/2+1;
    
    if da_sig>0
        dt=lags_xc(2)-lags_xc(1);
        a=1/(2*da_sig); %heigth of distribution
        sig_shift=2*da_sig/dt;
        s=-2*da_sig:dt:2*da_sig;
        f=zeros(size(s));
        f(s<=0)=s(s<=0)*a/(2*da_sig)+a;
        f(s>0)=-s(s>0)*a/(2*da_sig)+a;
        utent=u0-sig_shift:u0+sig_shift;
    else
        utent=u0;
        f=1;
    end
    
    for it=1:length(pre)
        C_delta(post(it),post(it),utent)=expected_count(pre(it))*f+squeeze(C_delta(post(it),post(it),utent))';
    end  
    
    %running C_delta through the expansion
    C_correction=C_delta;
%     for k=1:n
%         i_j_array=allVL1(2,k); %all possible pair of integers (i,j) summing to k
%         for k_runner=1:size(i_j_array,1) %running over all combos of i and j summin to k
%             %extracting current values of i and j
%             i=i_j_array(k_runner,1);
%             j=i_j_array(k_runner,2);
%             %extracting shifted C_hat
%             t_shift=(i-j)*da; %current time-shift
%             C_temp=dshift(C_hat,t_shift,dt); %shifted copy of C_hat
%             for u=1:size(C,3);
%                 C_correction(:,:,u)=C_correction(:,:,u)+J^i*squeeze(C_temp(:,:,u))*J'^j;
%             end
%         end
%     end

    %adding correction to C estimate
    C=C+C_correction;
end

end %function C_hat2C

function [A]=dshift(B,t_shift,dt)
    d_shift=t_shift/dt;
    index=1:size(B,3);
    index=index+d_shift;
    index(index>length(index))=length(index);
    index(index<1)=1;
    A=B(:,:,index);
end
    