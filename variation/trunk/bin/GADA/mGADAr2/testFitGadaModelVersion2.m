%TESTFITGADAMODEL Test the GADA-JRN and GADA-SMN algorithms
clear all;
disp('# GADA - Genome Alteration Detection Analysis');
disp('# Copyright (C) 2008,2009 Childrens Hospital Los Angeles');
disp('# Author: Roger Pique-Regi rpique@gmail.com');

%% Read data form text file... 
% It is better to write a small script to dump line by line to a binary
% file, instead of loading everything to memory and then dumping to the
% binary file

Ymat = dlmread('smallexample.txt', '\t');
Ymat = Ymat';
[N,M]= size(Ymat);

%% Placing input data in a binary file for memmory mapping
% The first step to use GADA-JRN is to create a binary file in which 
% the data is dumped. This file then is mapped in memmory, so we get more
% memory efficiency if very large files are used.

fd=fopen('mytemp.tmp','wb');
fwrite(fd,Ymat','double');
fclose(fd);
Ymat2=memmapfile('mytemp.tmp','format',{'double',[M N], 'x'});

%% Parameter settings
% See more information on these parameters with "help fitGadaModel"
% The most important are the following:
aAlpha=0.5;   %Increase this value and the algorithm will be faster miss some breakpoints (default 0.5)
TBackElim=4;  %Increase this value for less breakpoints ...
MinSegLen=3;  %Minumum number of probes for calling a segment.
verbose=2;    %Show debugging messages.

CommonParameters={...
    'TBackElim',TBackElim,'MinSegLen',MinSegLen,...
    'verbose',verbose};

%  isRemovingOK=0;
%  CommonParameters={...
%      'TBackElim',TBackElim,'MinSegLen',MinSegLen,...
%      'isRemovingOK',isRemovingOK,...
%      'AdjustSampleSigma2',0,'SampleSigma2',sigma2vec...
%      'verbose',verbose};

%% oldGADA 
[OldGADA.Iext,OldGADA.Wext,OldGADA.par]=fitGadaModel(Ymat2,...
    'aAlpha',aAlpha,'RefAdjust',0,'RhoAdjust',0,...
    CommonParameters{:});

%% GADA_JRN
% [GADA_JRN.Iext, GADA_JRN.Wext, GADA_JRN.par]=fitGadaModel(Ymat2,...
%     'aAlpha',aAlpha,'RefAdjust',1,'RhoAdjust',0,'TBackElim',3,'aAlpha',0.2,'MinSegLen',0,...
%     CommonParameters{:});

%% GADA_JRN
[GADA_JRN.Iext, GADA_JRN.Wext, GADA_JRN.par]=fitGadaModel(Ymat2,...
    'aAlpha',aAlpha,'RefAdjust',1,'RhoAdjust',1,'TBackElim',3,'aAlpha',0.2,'MinSegLen',0,...
    CommonParameters{:});

%% GADA_SMN
[GADA_SMN.Iext,GADA_SMN.Wext,GADA_SMN.par]=fitGadaModel(Ymat2,...
    'aAlpha',aAlpha,'RefAdjust',2,'RhoAdjust',0,...
    CommonParameters{:});

%% We don't need Ymat2 more, 
clear Ymat2

%% Final results Plotting

figure(1);clf;

% Original (true seg)
% subplot(2,3,1);
% imagesc(Xmat,[-2.0 2.0]);
% xlabel('Probe $m$','Interpreter','latex','fontsize',10);
% ylabel('Sample $n$','Interpreter','latex','fontsize',10);
% title('Underlying CNA component $x$','Interpreter','latex','fontsize',11);
% text(-M*0.15,-1,'A','fontsize',14);

% Observed intensity (Ymat)
%figure(5);
subplot(2,3,2);
imagesc(Ymat);
%imagesc(Ymat);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('Observed array data $y_{mn} = x_{mn} + \rho_n r_m + \epsilon_{mn}$','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'B','fontsize',14);

%% GADA_SMN
subplot(2,3,4);
YmatN=zeros(N,M);
for n=1:N
    YmatN(n,:)=Ymat(n,:)-GADA_SMN.par.Ref;
end
imagesc(YmatN);
%imagesc(YmatN);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('Median normalized data: $\tilde y = y -$ median $(y)$','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'D','fontsize',14);

%% GADA_JRN
subplot(2,3,6);
Xrec=zeros(N,M);
for n=1:N
    Xrec(n,:)=ReconstructFromIextWext(GADA_JRN.Iext{n},GADA_JRN.Wext{n});
end

imagesc(Xrec);
%imagesc(Xrec);

xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('GADA with joint reference normalization (GADA-JRN)','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'F','fontsize',14);

%% 
subplot(2,3,3);
Xrec2=zeros(N,M);
for n=1:N
    Xrec2(n,:)=ReconstructFromIextWext(OldGADA.Iext{n},OldGADA.Wext{n});
end
imagesc(Xrec2);
%imagesc(Xrec2);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('GADA without normalization','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'C','fontsize',14);

%% 
subplot(2,3,5);
Xrec3=zeros(N,M);
for n=1:N
    Xrec3(n,:)=ReconstructFromIextWext(GADA_SMN.Iext{n},GADA_SMN.Wext{n});
end
imagesc(Xrec3);
%imagesc(Xrec3);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('GADA with separate median normalization (GADA-SMN)','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'E','fontsize',14);

%%

delete *.tmp

