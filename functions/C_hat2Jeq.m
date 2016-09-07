%function that returns the estimate J equilibrium matrix
%as well as the associated stationary cross-corr, for a C_hat and some
%parameters p.

function [Jeq,C]=C_hat2Jeq(expansion,W_rule,J_init,lags_xc,p)

%internal parameters
max_it=10;
tol=1e-7;
dt=lags_xc(2)-lags_xc(1);

%range of synaptic strength for D function
j_range=p.w_min:0.0001:p.w_max;

J1=J_init;
J2=10*ones(3);
Jdiff=norm(J1-J2);
it=1;

%Iteration to equilibrium
while it<=max_it && Jdiff>tol
    it=it+1;
    %extracting self-consistent xcorrs
    J=J1*p.net_connect_prob*p.N/3;
    C=expansion(J);    
    %computing DJ
    for pre=1:3
        for post=1:3
              D=zeros(1,length(j_range));
              y=squeeze(C(post,pre,:));
              for j=1:length(j_range)
                D(j)=sum(y'.*W_rule(lags_xc-p.dendrite_del+p.net_del,j_range(j)))*dt;
              end
              ind=find(D<=0,1,'first');
%               figure;plot(D); pause(0.1)
              if isempty(ind)
                  J2(post,pre)=p.w_min;
              else
                  J2(post,pre)=(j_range(ind)+j_range(ind-1))/2;
              end
        end
    end
    Jdiff=norm(J2-J1);
    J1=J2;
end
display(['Found equilibrium J in ' num2str(it) ' iterations with ' num2str(Jdiff,10) ' error'])
Jeq=J2;
end

