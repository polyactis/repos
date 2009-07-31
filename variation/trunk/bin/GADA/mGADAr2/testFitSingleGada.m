%TESTFITSINGLEGADA Test the single sample GADA algorithm
disp('# GADA - Genome Alteration Detection Analysis');
disp('# Copyright (C) 2008,2009 Childrens Hospital Los Angeles');
disp('# Author: Roger Pique-Regi rpique@gmail.com');

L2=200;
sigma2=1;
L=100;
N=L*L2;

freqchirp=0.010;

origseq=sign(chirp(1:N,freqchirp,N/2,0.00));
origseq(L2/2:L2:end)=5  %Creat single probe outliers
origseq=origseq-mean(origseq);

dispos=find(diff(origseq)~=0);

testseq=origseq+randn(1,N)*sqrt(sigma2);
testseq=testseq-mean(testseq);

%%
%a=1.2
a=0.2
T=4;

[Iext,Wext]=fitSingleGada(testseq,sigma2,a,T,0);
tr1=ReconstructFromIextWext(Iext,Wext);

%[Iext2,Wext2]=fitSingleGada(testseq,sigma2,a,T,3);
[Iext2,Wext2]=BackwardElimination(Iext,Wext,sigma2,T,3);
tr2=ReconstructFromIextWext(Iext2,Wext2);

%% 
figure(2);
plot(testseq,'r.');
hold on;
plot(origseq,'g');
plot(tr1,'m','LineWidth',2);
plot(tr2,'b','LineWidth',2);

hold off;
