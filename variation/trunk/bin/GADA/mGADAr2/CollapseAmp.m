function [OutAmp,State]=CollapseAmp(SegAmp,SegLen,BaseAmp,sigma2,T)

% CollapseAmp -- Collapses the amplitudes of the non significantly altered segments into a base level. 
% [OutAmp,State]=CollapseAmp(SegAmp,SegLen,BaseAmp,sigma2,T);
% Matlab implementation, 
% See also: MexCollapseAmpTtest

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


OutAmp=SegAmp;
SegLen=double(SegLen);

OutAmp(CompTscore(BaseAmp,SegAmp,SegLen,sigma2)<T)=BaseAmp;
if nargout>1
    State=zeros(size(OutAmp)); %Neutral
    State(OutAmp>BaseAmp)=1; %Gain
    State(OutAmp<BaseAmp)=-1; %Loss
end

function T=CompTscore(m1,m2,n2,sigma2)

T=abs(m1-m2)./sqrt(sigma2./n2);


