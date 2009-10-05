function [Iext,Wext,par]=fitGadaModel(Ymat,varargin)

% fitGadaModel - Fits the original GADA model and its extensions to a large set of arrays. 
%
%   [Iext,Wext,par]=fitGadaModel(Ymat,...)
%
%   This function fits a GADA (Genome Alteration Detection Analysis) model 
% proposed in [1] and the extended model for joint normalization proposed 
% in [2]. This implementation uses memory mapped files to achieve analyzing 
% a large collection of samples with high density arrays with limited 
% physical memory. For the old single sample algorithm [1], the MexSBL_BE 
% is faster and does not require memory mapped files. 
% 
% The  GADA model for the previous version [1] and the new GADA-JRN [2] 
% (GADA with Joint Reference Normalization) can be written toghether as:
%       y_mn = x_mn + Rho_n * r_m + e_mn
%  where:
%    - y_mn are the observed array probe hybridization intensities, ordered
%      by chromosomal position. 
%    - x_mn is the underlying copy number effect on the array intensity. 
%    - r_m specifies a probe specific hybridization effect 'Ref' and Rho_n is 
%      sample specific scale factor for r_m. The old GADA model [1] would
%      assume Rho_n=0 and r_m=0. The new GADA-JRN [2] can estimate these
%      parameters jointly with the copy number. 
%    - e_mn is the noise assumed to be additive white gaussian N(0,sigma2_n). 
%      The variance of the noise sigma2_n can also be estimated inside
%      GADA. 
%  The copy number component x_mn is represented using a Piece-Wise
%  Constant (PWC) vector representation, for each sample 'n':
%             x^n = Fw^n
%  where F is a matrix with normalized step vectors and w^n is a SPARSE
%  vector of weights with a one to one correspondence to each copy number
%  breakpoint. The w^n are modeled using a Sparse Bayesian Learning
%  hierarchical prior:
%             w_mn|alpha_mn are i.i.d. N(0,1/alpha_mn)
%             alpha_mn are i.i.d Gamma(aAlpha,bAlpha)
%
% Required input parameters:
%   Ymat <- Observed, MxN matrix (N = no samples, M = no of probes). 
% Assumed to be a memmaparray object. testFitGadaModel.m shows how these 
% objects can be prepared and more information can be found on memmapfile
% help. 
% 
% The algorithm returns a set of breakpoint positions Iext and the
% magnitude of the breakpoints Wext vector for each of the samples and are
% gathered toghether using a Matlab cell array.
% 
%    All the other input arguments are optional and given as a "'Property', 
% value" pair to the algorithm; e.g. fitGadaModel(...,'Property',value,...)
%
% 'SampleSigma2' <- Variance of the noise present on each sample (vector of size
%   N). Note that the algorithm can also estimate these parameters with
%   'AdjustSampleSigma2=1'. If this estimation is activated then Sigma2 is taken
%   as initial point on the EM algorithm otherwise 'Sigma2',ones(N,1) is
%   used by default. 
%
% 'aAlpha' <-- Prior distribution hyperparameters, ('aAlpha',0.5) default. 
% 'bAlpha' <-- Prior distribution hyperparameters, ('bAlpha',0.0) default. 
%    bAlpha does not need to be changed, a higher aAlpha implies a more sparse 
%    prior and less breakpoints will be reported, leading to a lower FDR. 
% 
% 'MaxAlpha','isRemovingOK' <- Maximum value for the Alpha parameters before
%    we consider we can remove it from the model if isRemovingOK=1.
%    It usually gives a speed-up as more breakpoints are removed. The default
%    values are 'MaxAlpha',1E7 and 'isRemovingOK',1. 
%
% 'MaxIt' <- Maximum number of EM iterations before we stop waiting for the
%    algorithm to converge (1000). 
% 'Tol' <- Minimum relative change in the objective function before we
%    consider the algorithm has converged. 
%
% 'verbose' <- show verbose messages, 0 no messages, 1 or higher more
%    messages are returned. 
%
% Parameters for GADA-JRN [2]: 
% 'Ref' <- Reference intensity 'r_m' or Probe specific hybridization bias. 
%    (default r_m=0, Ref=zeros(1,M), i.e. and no reference is removed).
% 'RefAdjust' <- Sets how 'Ref' is ajusted/corrected:
%    If 'RefAdjust=0' 'Ref' is fixed, the 'Ref' provided/default is used.  
%    If 'RefAdjust=1' 'Ref' is re-estimated inside the EM loop (GADA-SMN)
%                   starting with the provided/default value for 'Ref'.
%    If 'RefAdjust=2' 'Ref' is fixed, calulated as the median across samples
%                   before the aglorithm starts (GADA-SMN)          
% 'Rho' <- Scale parameter for the Ref artifact (default Rho=ones(1,N)). 
% 'RhoAdjust' <- 0 means that initial 'Rho' is used all the time, 1 means
%    that 'Rho' is ajusted in the EM algorithm. 
%  
% Backward elimination parameters [1,2]:
% 'TBackElim' <- Critical value T of the Backward elimination. All the
%    remaining breakpoints should have a score larger than T. The score
%    has a null distribution in case of Gaussian noise and it is basically
%    the difference of the averages between to neighboring segments divided
%    by the standard error of that difference (taking into account the
%    noise variance and the size in probes of the neighboring segments).
%    Default TBackElim=4
% 'MinSegLen' <- Minimum size of a segment in probes (MinSegLen=3).
%    Segments with less probes will be considered outliers and removed.
%
% References:
% [1] Pique-Regi R, Monso-Varona J,Ortega A, Seeger RC, Triche TJ,
%   Asgharzadeh S: "Sparse representation and Bayesian detection of the
%   genome copy number alterations from microarray data", Bioinformatics ,
%   Feb 2008
% [2] Pique-Regi R, Ortega A, Asgharzadeh S: "Joint estimation of copy
%   number variation and reference intensities on multiple DNA arrays using
%   GADA", Bioinformatics , Feb 2009
%
% See also:
%   ReconstructFromIextWext.m IextWextToSegAmp.m fitSingleGADA.m

% %   This File is part of GADA
% % 
% %   GADA v2.0 Genome Alteration Detection Algorithm 
% %   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
% % 
% %   GADA is free software: you can redistribute it and/or modify
% %   it under the terms of the GNU General Public License as published by
% %   the Free Software Foundation, either version 3 of the License, or
% %   any later version.
% % 
% %   GADA is distributed in the hope that it will be useful,
% %   but WITHOUT ANY WARRANTY; without even the implied warranty of
% %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% %   GNU General Public License for more details.
% % 
% %   You should have received a copy of the GNU General Public License
% %   along with GADA.  If not, see <http://www.gnu.org/licenses/>.
% % 
% %   Authors: 
% %   Roger Pique-Regi    rpique@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parsing the input parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required parameters Ymat
[M,N]=size(Ymat.Data.x);

% Defining property parameters and default values.
par.verbose=1;
par.Tol=1E-8;
par.MaxIt=10000;
par.isRemovingOK=1;
par.MaxAlpha=1E7;
par.bAlpha=1E-20;
par.aAlpha=0.5;
par.AdjustSampleSigma2=1;
par.SampleSigma2=[];
par.RefAdjust=0;
par.Ref=[]; %InitialRef...
par.PreInitRefActivated=1; %Preinitialization of r
par.PreInitRefWinSize=100;
par.RhoAdjust=0; %Adjustment of different waveforms
par.Rho=[]; %Initial Rho
par.TBackElim=4;
par.MinSegLen=3;


%Parsing parameters
par=parsepvpairsparstruct(par,varargin{:});

Kglobal=M-1; %Number of breakpoints basis elements

if(par.verbose>0)
    fprintf('Size of the input data: \n')
    fprintf('   M=%d, N=%d, K=%d\n',M,N,Kglobal);
    fprintf('Running %s settings:\n',mfilename)
    disp(par)
    fprintf('\n');
end

if(par.RhoAdjust==1 & par.RefAdjust~=1)
    error('RhoAdjust=1 requires RefAdjust=1');
end
if(par.bAlpha>par.aAlpha/100)
    warning(mfilename,'bAlpha should be 0 or a very small value');
end
if(par.aAlpha<0.1)
    warning(mfilename,'EM may take a much longer time to converge for small values of aAlpha, typically aAlpha>0.1');
end
if(par.aAlpha>2)
    warning(mfilename,'Sparse Bayesian Learning hyperprior may not work well for large values of aAlpha');
end
assert(par.aAlpha>=0);
assert(par.bAlpha>=0);
assert(par.TBackElim>=0);
assert(par.MinSegLen>=0);
assert(par.AdjustSampleSigma2==1 || par.AdjustSampleSigma2==0);
assert(par.RhoAdjust==1 || par.RhoAdjust==0);
assert(par.RefAdjust==0 || par.RefAdjust==1 || par.RefAdjust==2);
assert(par.MaxAlpha>1000);
assert(par.MaxIt>=0);
assert(par.Tol>0 && par.Tol<0.1);
assert(all(par.SampleSigma2>0));
assert(all(par.Rho>0));


%%%%%%%%%%%%%%%%%%%%%%
%% Initializations  %%
%%%%%%%%%%%%%%%%%%%%%%


%Variance estimation on sample level
if isempty(par.SampleSigma2)
    %TODO estimate Sigma2
    par.SampleSigma2=ones(1,N);
else
    par.SampleSigma2=par.SampleSigma2;
end
%aBeta=M/20; %Trust in the initial estimation...
aBeta=100/2; %Trust in the initial estimation... it should not be important if there is a large number of probes. 
bBeta=par.SampleSigma2.*aBeta;

%% Initialize memmapfile for the large variables in the model...
%We will start with stuff that need to be kept between EM iterations.
internalmapfile=['temp',num2str(round(now*100000)),'.tmp'];
verbose(par.verbose,' Preparing memmapfile %s\n',internalmapfile);
if exist(internalmapfile,'file')
    delete(internalmapfile);
end
fd=fopen(internalmapfile,'wb');
%Adding Iext=int32(0:M) N times
for n=1:N 
    fwrite(fd,int32(0:M),'int32'); 
end
%Adding Alphas=zeros(K,N)
for n=1:N 
    fwrite(fd,zeros(Kglobal,1)+1E-20,'double'); 
end
%Adding    W0=zeros(N,M-1);
for n=1:N 
    fwrite(fd,zeros(Kglobal,1),'double'); 
end
fclose(fd);
%Tscores=zeros(size(Alpha)); 

m=memmapfile(internalmapfile,'format',...
    {'int32' [Kglobal+2 N] 'Iext';...
     'double' [Kglobal N] 'Alpha';...
     'double' [Kglobal N] 'W0'},'writable',true);

%% Initialize h
verbose(par.verbose,'Initalizations \n');
[h0,h1]=CompH(M); %PWC dependent....

Kvec=ones(1,N)*Kglobal;

%Ref initialization
if isempty(par.Ref)
    par.Ref=zeros(1,M); 
else
    assert(length(par.Ref)==M);
end

if par.RefAdjust==2
    verbose(par.verbose,'Estimating probe hybridization with median (GADA-SMN)');    
    BlockSize=round(32E6/N);
    mEnd=1;
    while mEnd<M
        mIni=mEnd;
        mEnd=min(M,mEnd+BlockSize);
        %Block=Ymat.Data.x(mIni:mEnd,:)
        par.Ref(mIni:mEnd)=median(Ymat.Data.x(mIni:mEnd,:),2);
    end
end
    
if par.RefAdjust==1 & par.PreInitRefActivated==1
    verbose(par.verbose,'Pre-Initial Ref estimation with window method\n');        
    Ref=zeros(1,M);
    mStart=1;
    for mEnd=1:par.PreInitRefWinSize:M
%        progress(,M);
        WinStart=max(1,mStart-par.PreInitRefWinSize);
        WinEnd=max(M,mEnd-par.PreInitRefWinSize);
        
        WindowMedianPerSample=median(Ymat.Data.x(mStart:mEnd,:)); 
        WindowMedian=median(WindowMedianPerSample);
        verbose(par.verbose-1,'Probes: %d --%d    Median: %g \n',mStart,mEnd,WindowMedian);                
        Ref(mStart:mEnd)=WindowMedian;
        mStart=mEnd+1;
    end
    fprintf('\n');
end

%Rho adjustment initialization
if isempty(par.Rho)
    par.Rho=ones(1,N);
else
    assert(length(par.Rho)==N);
end

% Averages, 
yavg=zeros(1,N);
for n=1:N
    yavg(n)=mean(Ymat.Data.x(:,n)); %% Columnwise avg.
end


% Initialize W_0, Z matrix ...
%    Zmat=zeros(N,M-1); %Zmat is necessary only for base removal,
for n=1:N %This is only necessary once if RefAdjust==0
    %m.Data.W0(:,n)=CompFdualXb(Ymat.Data.x(:,n)'-par.Ref-yavg(n)); %Remove the avg and ref for each sample
    %        Zmat(n,:)=TridSolve(h1,h0,h1,W0rMat(n,:));    %Compute the z vectors for each sample
    myZ=MexCompZ_PWCnorm(Ymat.Data.x(:,n)'-par.Rho(n)*par.Ref-yavg(n));  %PWC dep
    m.Data.W0(1:Kvec(n),n)=TriSymGaxpy(h0,h1,myZ(m.Data.Iext(2:(Kvec(n)+2-1),n)));
    
end 

par.LogLikEvolution=zeros(1,par.MaxIt);
LogLikSampleVec=zeros(1,N);

%%%%%%%%%%%%%%%%%%%%%
%% First round SBL %%
%%%%%%%%%%%%%%%%%%%%%
%[WrMat,Iext,Alpha,SigmaWr]=itera_Multi_SBL_PWC(Zmat,par.SampleSigma2,aAlpha,par.bAlpha,Iext,Alpha,par.MaxAlpha,isRemovingOK,MaxIt,Tol,par.verbose);

par.SingleBrkRemovals=0; 
par.MultiBrkRemovals=0;

verbose(par.verbose,'Running EM loop \n');
for i=1:par.MaxIt
    
    if par.verbose>1
        %progress(i);
        fprintf('====================================================\n');
        fprintf('Iteration Number: %d\n',i);
    end
    
    if par.RefAdjust==1
        rnew=zeros(1,M);
        if par.RhoAdjust==1
            DenR=sum(par.Ref.^2); 
        end
    end

    for n=1:N
        K=Kvec(n);
        Alpha=m.Data.Alpha(1:K,n);
        
        verbose(par.verbose-2,' Processing Sample %d\n',n);                        
        %%%%%%%%%%%%%%%%%%%%%%
        %%%  Subselection   %%
        %%%%%%%%%%%%%%%%%%%%%%
        if (par.isRemovingOK)
            %            sel=find(sum((m.Data.Alpha(1:K,:))<par.MaxAlpha,2)>0); %Find global breakpoints.
            sel=find(Alpha<par.MaxAlpha);
            Iext=m.Data.Iext(1:K+2,n);
            if length(sel)<K
                if (K-length(sel))==1
                    par.SingleBrkRemovals=par.SingleBrkRemovals+1;
                else
                    par.MultiBrkRemovals=par.MultiBrkRemovals+1;                   
                end
                verbose(par.verbose-2,'  Removing %d breakpoints of %d on sample %d\n',K-length(sel),K,n);
                Iext = Iext([1 sel'+1 K+2]);
                %%%% we need to do it all the time now since each sample
                %%%% can have different alphas.
                %[h0,h1]=CompHs1(double(Iext(2:end-1)),M);  % PWC dep
                K=length(Iext)-2;
                Alpha=Alpha(sel); %m.Alphas will be done in the Mstep
                m.Data.Iext(1:K+2,n)=Iext;
                %W0rMat(n,:)=TriSymGaxpy(h0,h1,Zmat(n,Iext(2:end-1)));
                Kvec(n)=K;
            end
            if K>0
                [h0,h1]=CompHs1(double(Iext(2:end-1))',M);  % PWC dep
            else %No breakpoints left on this sample...                
                if par.RefAdjust==1
                    verbose(par.verbose-2,'  Partial M-Step for reference \n');                    
                    if par.RhoAdjust==1
                        %new Rho
                        if(abs(DenR)>1E-3)
                            par.Rho(n)=sum(Ymat.Data.x(:,n)'.*par.Ref)/DenR;
                        end
                        rnew=rnew+par.Rho(n)*(Ymat.Data.x(:,n)'-yavg(n)); %gplus...
                    else                        
                        rnew=rnew+(Ymat.Data.x(:,n)'-yavg(n)); %gplus...
                    end
                    %%% Remember to divide par.Ref by N !!!
                end      
                %%%  Reference remove   %%
                if par.RhoAdjust 
                    UnrefY=Ymat.Data.x(:,n)'-par.Rho(n)*par.Ref-yavg(n);
                else                                            
                    UnrefY=Ymat.Data.x(:,n)'-par.Ref-yavg(n);
                end
                RSS=sum((UnrefY).^2);
                if par.AdjustSampleSigma2==1
                    verbose(par.verbose-2,'  M-Step for variance on sample %d,',n);
                    verbose(par.verbose-2,' Old = %g', par.SampleSigma2(n));          
                    par.SampleSigma2(n)=(RSS+2*bBeta(n))/(M+2*aBeta);
                    verbose(par.verbose-2,' New Sigma2 = %g\n', par.SampleSigma2(n));
                end
                verbose(par.verbose-2,'  Partial LogLik for sample %d: LogLik = ',n);
                LogDetSigmaY=(M-K)*log(par.SampleSigma2(n));
                LogLikSample=-1/2*LogDetSigmaY-1/2*(RSS/par.SampleSigma2(n));%-M/2*log(2*pi);
%                LogLikSample=LogLikSample+(par.aAlpha)*sum(log(Alpha))-par.bAlpha*sum(Alpha);
                LogLikSample=LogLikSample-(aBeta)*log(par.SampleSigma2(n))+bBeta(n)*(1/par.SampleSigma2(n));
                LogLikSample=LogLikSample+(M-1-K)*((par.aAlpha)*log(par.MaxAlpha) - par.bAlpha*par.MaxAlpha);
                verbose(par.verbose-2,'LogLik Sample %d is %g Change %g \n',n,LogLikSample,LogLikSample-LogLikSampleVec(n));
                LogLikSampleVec(n)=LogLikSample;
                par.LogLikEvolution(i)=par.LogLikEvolution(i)+LogLikSample; 
                continue;
            end %No breakpoints left on this sample...
        else
            %If we don't remove bases we should still control overflow of the alphas. 
            sel=find(Alpha>par.MaxAlpha);
            Alpha(sel)=par.MaxAlpha*100; %I don't want to let it go to a very big value to overflow.
            verbose(par.verbose-2,'  Removing %d breakpoints of %d on sample %d\n',length(sel),K,n);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  Reference remove   %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        verbose(par.verbose-2,' Substracting reference for sample %d \n',n);
        if par.RhoAdjust
            %new Rho
            if(abs(DenR)>1E-3)
                par.Rho(n)=sum(Ymat.Data.x(:,n)'.*par.Ref)/DenR;
            end
            %Corrected data.
            UnrefY=Ymat.Data.x(:,n)'-par.Rho(n)*par.Ref-yavg(n);
        else
            UnrefY=Ymat.Data.x(:,n)'-par.Ref-yavg(n);
        end

        if par.isRemovingOK
            myZ=MexCompZ_PWCnorm(UnrefY);
            m.Data.W0(1:K,n)=TriSymGaxpy(h0,h1,myZ(Iext(2:(K+2-1))));
        else
            m.Data.W0(:,n)=CompFdualXb(UnrefY);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%       E STEP       %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        verbose(par.verbose-2,'  E-Step on sample %d \n',n);
        [SigmaWn,MeanWn]=comp_Estep_SBL(m.Data.W0(1:K,n)',h0,h1,Alpha',par.SampleSigma2(n));
       
        if (par.isRemovingOK)
            ReconY=ReconstructFromIextWext(Iext,[0 MeanWn]);
        else
            ReconY=ReconstructFromIextWext(int32(0:M),[0 MeanWn]);            
        end
        RSS=sum((UnrefY-ReconY).^2); % for the noise and for the L

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%       M STEP       %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%% Noise Update %%%%%%%%%%%%%%%%%
        %M-step, for the noise        
        if par.AdjustSampleSigma2==1
            verbose(par.verbose-2,'  M-Step for variance on sample %d,',n);
            verbose(par.verbose-2,' Old = %g', par.SampleSigma2(n));
            
            SumGamma=sum(1-Alpha.*SigmaWn');
            %RSS=sum((UnrefY-ReconY).^2);
            par.SampleSigma2(n)=(RSS+par.SampleSigma2(n)*SumGamma+2*bBeta(n))/(M+2*aBeta);
            verbose(par.verbose-2,' New Sigma2 = %g\n', par.SampleSigma2(n));
            
            if(isnan(par.SampleSigma2(n)))
                error('NaN occurred');
            end
            
        end        

        %%%%%%%%%  Alpha Update %%%%%%%%%%%%%%%%
        verbose(par.verbose-2,'  M-Step for alpha on sample %d:',n);
        Alpha=(1+2*par.aAlpha)./(MeanWn.^2+SigmaWn+2*par.bAlpha);

        verbose(par.verbose-2,' Change max (1/alpha-1/alpha*) = %g \n',max(abs(1./(m.Data.Alpha(1:K,n))-1./Alpha'))); 
        m.Data.Alpha(1:K,n)=Alpha;       
        

        %%%%%%%%% Reference update %%%%%%%%%%%%%%%
        if par.RefAdjust==1
            verbose(par.verbose-2,'  Partial M-Step for reference \n');
            %par.Ref=par.Ref+(Ymat(n,:)-yavg(n)-ReconstructFromIextWext(Iext,[0 WrMat(n,:)]));
            if par.RhoAdjust==1
                rnew=rnew+par.Rho(n)*(Ymat.Data.x(:,n)'-yavg(n)-ReconY); 
            else
                rnew=rnew+(Ymat.Data.x(:,n)'-yavg(n)-ReconY); 
            end
            %%% Remember to divide par.Ref by N !
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%       LogLik compute          %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        verbose(par.verbose-2,'  Partial LogLik for sample %d \n ',n);
        [xF,Fu,Fd,Fl] = MexForwardTridSolve(h1,h0,h1,h0);
        LogDetH=sum(log(Fd));
        t1=par.SampleSigma2(n)*h1;
        t0=par.SampleSigma2(n)*h0+1./Alpha; 
        [xF,Fu,Fd,Fl] = MexForwardTridSolve(t1,t0,t1,t0);
        LogDetT2=sum(log(Fd));
        LogDetSigmaY=LogDetT2-LogDetH+(M-K)*log(par.SampleSigma2(n));
        ASS=sum(Alpha.*MeanWn.*MeanWn);
        LogLikSample=-1/2*LogDetSigmaY-1/2*(RSS/par.SampleSigma2(n)+ASS);%-M/2*log(2*pi);
        LogLikSample=LogLikSample+(par.aAlpha)*sum(log(Alpha))-par.bAlpha*sum(Alpha);
        LogLikSample=LogLikSample-(aBeta)*log(par.SampleSigma2(n))+bBeta(n)*(1/par.SampleSigma2(n));
        LogLikSample=LogLikSample+(M-1-K)*((par.aAlpha)*log(par.MaxAlpha) - par.bAlpha*par.MaxAlpha); %To account already removed stufff.

        if i>1
            verbose(par.verbose-2,'LogLik Sample %d is %g Change %g \n',n,LogLikSample,LogLikSample-LogLikSampleVec(n));
%            if LogLikSample-LogLikSampleVec(n) < 0           
%                warning('Log lik decrease');
%            end
%       NOTE: When we use isRemovingOK sometimes we can have a false
%       negative increase of the likelihood due to how we do the
%       calculations, but this is not a problem. 
        else
            verbose(par.verbose-2,'%g \n',LogLikSample);
        end
        if(isnan(LogLikSample))
            %error('');
            return
        end
        LogLikSampleVec(n)=LogLikSample;
        
        par.LogLikEvolution(i)=par.LogLikEvolution(i)+LogLikSample; %gplus

    end %% Loop along samples 

    if par.RefAdjust==1
        verbose(par.verbose-1,' Final M-Step for reference: Sum(abs(Ref-Ref*)) change=%g  \n',mean(abs(par.Ref-rnew/N)));                
        if par.RhoAdjust==1
            verbose(par.verbose-1,' Final M-Step for Rho:\n');                
            par.Ref=rnew/sum(par.Rho.^2);
        else
            par.Ref=rnew/N;                    
        end
    end
    
    if (par.isRemovingOK) && (par.verbose>1)
        fprintf('Breakpoints:');
        fprintf(' %d',Kvec);                        
        fprintf('\n');        
    end
    if (par.AdjustSampleSigma2) && (par.verbose>1)
        fprintf('Sigma2:');
        fprintf('%g ',par.SampleSigma2);                        
        fprintf('\n');        
    end

    
    verbose(par.verbose-1,' Total Log-Lik iteration %d is %g\n',i,par.LogLikEvolution(i))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Determine convergence   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if i>1
        verbose(par.verbose-1,' Increase versus previous iteration %g\n',par.LogLikEvolution(i)-par.LogLikEvolution(i-1));
        if(abs((par.LogLikEvolution(i)-par.LogLikEvolution(i-1))/par.LogLikEvolution(i))<par.Tol)
            verbose(par.verbose,'====================================================\n');
            verbose(par.verbose,'\n yes! the model converged in %d itarations!\n',i);
            break;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fancy verbose stuff, with graphical plots %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    if par.verbose>50
%          verbose_v1(Iext,Alpha,MaxAlpha,N,M,WrMat,SigmaWr,mynorm,i);
%        fancy_verbose(m,par);
%    end
    %%%%%%%%%%%%%%%%%%%%%
end %EM loop

par.LogLikEvolution=par.LogLikEvolution(1:i);

par.NumIt=i;

if i>=par.MaxIt
    verbose(par.verbose,'Converged??? Stopped after %d iterations with change %g\n',i,par.LogLikEvolution(i)-par.LogLikEvolution(i-1));
end


if (par.isRemovingOK) && (par.verbose>0)
    verbose(par.verbose,'Total Single Removals =%d, Multiple Removals =%d \n',par.SingleBrkRemovals,par.MultiBrkRemovals);
    fprintf('Breakpoints:');
    fprintf(' %d',Kvec);
    fprintf('\n');
end
if (par.AdjustSampleSigma2) && (par.verbose>0)
    fprintf('Sigma2:');
    fprintf('%g ',par.SampleSigma2);
    fprintf('\n');
end

%% Backward elimination and final breakpoint report.

Wext=cell(1,N);
Iext=cell(1,N);
%Alpha=cell(1,N);
%SigmaW=cell(1,N);

verbose(par.verbose,'Computing final weights\n');
for n=1:N
    K=Kvec(n);
    if par.isRemovingOK
        Iext{n}=m.Data.Iext(1:K+2,n);
        if (K>0)
            [h0,h1]=CompHs1(double(m.Data.Iext(2:K+1,n))',M);  % PWC dep
        end
    else
        Iext{n}=int32(0:M); %or 0:M;
    end

    verbose(par.verbose-1,'Retrieving EM fit sample %d',n);
    if par.isRemovingOK
        myAlpha=zeros(1,K); %% Re-project 
    else
        myAlpha=m.Data.Alpha(1:K,n)'; 
        myAlpha(myAlpha<1000)=0.0; %% Re-project
    end
    if (K>0)
        [SigmaWn,MeanWn]=comp_Estep_SBL(m.Data.W0(1:K,n)',h0,h1,myAlpha,par.SampleSigma2(n));
    else
        MeanWn=[];
    end

    Wext{n}=[yavg(n), MeanWn];
    verbose(par.verbose-1,'  Backward elimination, ');
    [Iext{n},Wext{n}]=MexBEwTandMinLen(Iext{n},Wext{n},par.SampleSigma2(n),par.TBackElim,par.MinSegLen); 
    Kvec(n)=length(Iext{n})-2;
    verbose(par.verbose-1,'breakpoints removed = %d\n',K-Kvec(n));

 %   Alpha{n}=[];%m.Data.Alpha(1:K,n);
 %   SigmaW{n}=[];%SigmaWn;

end
if (par.isRemovingOK) && (par.verbose>0)
    fprintf('Breakpoints:');
    fprintf(' %d',Kvec);
    fprintf('\n');
end

%Think how I would do it if I wanted to keep m for other things. 

clear m; 

delete(internalmapfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h0,h1]=CompH(dim)

% CompH - Computes H=inv(F'F) Tridiagonal symetric F unit norm.
% [h0,h1]=CompH(dim)
% h0 Main diagonal
% h1 Lower/Upper diagonal
% (Normalized case)

h0=(dim-1:-1:1).*(1:dim-1)/dim;
h1=-sqrt(h0(2:dim-1).*h0(1:dim-2));
h0=h0*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h0,h1]=CompHs1(Ir,M)

% CompHs2 Computes Hs=inv(Fs'Fs) Tridiagonal symetric
%  [h0,h1]=CompHs1(Ir,M)
% Ir  Breakpoint Position Set
% M  length of the observed vector
% h0 Main diagonal
% h1 Lower/Upper diagonal
%
% This for the normalized basis case!
% 
% Extend this to the Iext notation.

Ms=length(Ir);
a=sqrt(Ir(1:Ms-1)./(M-Ir(1:Ms-1)).*(M-Ir(2:Ms))./Ir(2:Ms));
h0=1./(1-a.*a);
h1=-a.*h0;
h0(Ms)=1;
h0(2:Ms)=a.*a.*h0(1:Ms-1)+h0(2:Ms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sigma2Wr,MeanWr]=comp_Estep_SBL(W0r,h0,h1,Alpha,Sigma2)

% comp_Estep_SBL Implements the E_step of the SBL algorithm
% It can be used either in the multiple sample or the single sample
% approach.

K=length(h0);

%%% [t0,tu,tl]=CompT(h0,h1,Alpha,Sigma2);
t0=((h0(1:K).*Alpha(1:K)))+1/Sigma2;
tu=(h1(1:K-1).*Alpha(2:K));
tl=(h1(1:K-1).*Alpha(1:K-1));

if K==1
    MeanWr=W0r./t0*(1/Sigma2);
    Sigma2Wr=h0/t0;
else
    
    %%[x,err]= TridSolve(tu,t0,tl,W0r);
    [xF,Fu,Fd,Fl] = MexForwardTridSolve(tu,t0,tl,W0r);
    [xB,Bu,Bd,Bl] = MexBackwardTridSolve(tu,t0,tl,W0r);   
    
    MeanWr=(1/Sigma2)*0.5*(xB+xF);
    
    if(sum(abs(xF-xB))>1E-8)
        warning('!!!WARNING!!! System may not be solved acurately !!!WARNING!!!');
    end
    %compute diagonal of $\Sigma$
    %[itu,it0,itl] = TridOfInvTrid(tu,t0,tl);
    [itl,it0,itu] = MexTridOfInvTrid(Fu,Fd,Bl);

    %Sigma2Wr = DiagOfTriXTri(itl,it0,itu,h1,h0,h1); %MexDiagOfTriXTri
    Sigma2Wr = MexDiagOfTriXTri(itl,it0,itu,h1,h0,h1); %MexDiagOfTriXTri

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function fancy_verbose(m,par)

