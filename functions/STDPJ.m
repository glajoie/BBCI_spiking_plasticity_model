%function that gives the STDP rule

function [y]=STDPJ(delta_t,j,w_max,gamma,varargin)

    % t_post is a scalar
    % t_pre can be a vector
    
    %ARGUMENT GOES : delta_t=t_pre-t_post

    if ~isempty(varargin)
        d=varargin{1};
        c_p=d.c_p;
        tau_p=d.tau_p;
        c_d=d.c_d;
        tau_d=d.tau_d;
    else
            %PArameters from Gilson IV
        %   c_p=15;
        %   tau_p=17; %(ms)
        %   c_d=10;
        %   tau_d=34; %(ms)

        % Parameters from Gilson V
        c_p=30;
        tau_p=8.5;
        c_d=20;
        tau_d=17;
        % c_d=30;
        % tau_d=8.5;
    end

  %trim delta_t so it's not too big
%   delta_t=delta_t(delta_t>=-250);
%   delta_t=delta_t(delta_t<=250);
  y=zeros(size(delta_t));
  
  y(delta_t==0)=0;
  y(delta_t<0)=c_p*exp(delta_t(delta_t<0)/tau_p).*(-delta_t(delta_t<0))/tau_p*(1-j/w_max)^gamma;
  y(delta_t>0)=-c_d*exp(-delta_t(delta_t>0)/tau_d).*delta_t(delta_t>0)/tau_d*(j/w_max)^gamma;
    