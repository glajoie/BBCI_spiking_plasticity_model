%function that shifts cross-corr matrices by a delay ammount

function [A]=dshift(B,t_shift,dt)
    d_shift=t_shift/dt;
    index=1:size(B,3);
    index=index+d_shift;
    index(index>length(index))=length(index);
    index(index<1)=1;
    A=B(:,:,index);
end