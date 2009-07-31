function [C12,U1,U2]=comparedisc(q1,q2,w)

%comparedisc -- Compares two discontinuity sets using Willenbrock measure
% [C12,U1,U2]=comparedisc(q1,q2,w)
% q1 Breakpoint set 1
% q2 Breakpoint set 2
% w  Maximum window length to validate a breakpoint
% C12 1x(w+1) Brekpoints that are present in both q1 and q2 at any w=0:w
% U1 1x(w+1) Breakpoints that are only present at q1
% U2 1x(w+1) Breakpoints that are only present at q2
% Example:
% [auxTP,auxNF,auxFP]=comparedisc(realdisc,obsdisc,w);

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


if nargin<3
    w=0;
end

C12=zeros(1,w+1);
U1=zeros(1,w+1);
U2=zeros(1,w+1);
q1=sort(q1);
q2=sort(q2);


cw=0;
while (cw<=w)    
    if(cw>0)
        C12(cw+1)=C12(cw);
    end
    ind1=zeros(1,length(q1));
    for i1=1:length(q1)
        l2=length(q2);
        for i2=1:l2
            if(abs(q1(i1)-q2(i2))<=cw)
                C12(cw+1)=C12(cw+1)+1;
                q2=q2([1:(i2-1) (i2+1):l2]);
                ind1(i1)=1;
                break
            end
        end
    end    
    q1=q1(find(ind1==0));
    
    U1(cw+1)=length(q1);
    U2(cw+1)=length(q2);
    cw=cw+1;
end

    