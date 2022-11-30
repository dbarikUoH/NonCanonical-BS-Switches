% Bifurcation analysis of feedback refulated systems using potential
% energy function (PEF) calculation.
% Search of parameters for canonical and noncanonical bistable switches
% Count of different types of switches
% Segragation of parameters for different types of switches 

% Date: November 25, 2022
% Debashis Barik, School of Chemistry, University of Hyderabad, India

% Model: nnMA
% Logic Gate: OR

clear
clc
tic

% Initial values of the counts
dni=0;              % Dual isola
nibi=0;             % Bistable-Isola
iibi=0;             % Bistable-Inverted isola
nim=0;              % Isola-Mushroom
mbs=0;              % Bistable-Mushroom    
iim=0;              % Inverted isola-Mushroom
mush=0;             % Mushroom
invmush=0;          % Inverted mushroom
iso=0;              % Isola
inviso=0;           % Inverted isola
bis=0;              % Bistable
dbs=0;              % Dual bistable
nidbs=0;            % Dual bistable-Isola
tbs=0;              % Triple bistable
iidbs=0;            % Dual bistable-Inverted isola

% Array for searched parameters for differen switches
para_BS_all=[];
para_BS_pck=[];
para_mush=[];
para_invmush=[];
para_iso=[];
para_inviso=[];
para_bis=[];
para_dbs=[];
para_dni=[];
para_nibi=[];
para_iibi=[];
para_nim=[];
para_mbs=[];
para_iim=[];
para_nidbs=[];
para_tbs=[];
para_iidbs=[];

% To store the index of the parameters
dnipar=[];
nibipar=[];
iibipar=[];
invmushpar=[];
mushpar=[];
nimpar=[];
mbspar=[];
iimpar=[];
isopar=[];
invisopar=[];
bispar=[];
dbspar=[];
bistability_par=[];
nidbspar=[];
tbspar=[];
iidbspar=[];

ic1=0;
ic2=0;

rng('shuffle')

% Loading the thresholds obtained from half-function rules
load jas_in_nnMA_OR.dat
load jab_in_nnMA_OR.dat
load jba_in_nnMA_OR.dat

% Sample size of parameter combinations
nr=50000;                    

size_jas=length(jas_in_nnMA_OR);
size_jab=length(jab_in_nnMA_OR);
size_jba=length(jba_in_nnMA_OR);

% Randomizing the threshold parameters
rnd_as=randperm(size_jas)';
rnd_bs=randperm(size_jas)';
rnd_ab=randperm(size_jab)';
rnd_ba=randperm(size_jba)';

% Selecting parameters from uniform distributions
ga0=1+(10-1).*rand(nr,1);
gas=1+(100-1).*rand(nr,1);
gab=1+(100-1).*rand(nr,1);
gb0=1+(10-1).*rand(nr,1);
gbs=1+(100-1).*rand(nr,1);
gba=1+(100-1).*rand(nr,1);
nas=randi([1,10],nr,1);
nbs=randi([1,10],nr,1);
nab=randi([1,10],nr,1);
nba=randi([1,10],nr,1);
gma=0.01+(0.1-0.01).*rand(nr,1);
gmb=0.01+(0.1-0.01).*rand(nr,1);
% Selecting threshold parameters from distributions of half-functional rule
Jas=jas_in_nnMA_OR(rnd_as(1:nr));
Jbs=jas_in_nnMA_OR(rnd_bs(1:nr));
Jab=jab_in_nnMA_OR(rnd_ab(1:nr));
Jba=jba_in_nnMA_OR(rnd_ba(1:nr));

% Numerical parameters for integration: 
        % The maximum value of signal (S) for 1-p bifurcation run: 1000
sigE1=500;          % For low-res run: no. of points in the S scan          
sigE2=4000;         % For high-res run: no. of points in the S scan
dsig1=2.0;          % For low-res run: dS (interval of S)   
dsig2=0.25;         % For high-res run: dS (interval of S)
dx1=1;              % For low-res run: dB for integration of force over B
dx2=1;              % For high-res run: dB for integration of force over B
xE=10000;           % Maximum value of B    
b1=[0:dx1:xE];      
b2=[0:dx2:xE];
sig1=[1:sigE1];
sig2=[1:sigE2];

% Arrays for storing the steady states in the low-res run
nmbr_peaks_ss=ones(1,sigE1)*NaN;
nmbr_peaks_us=ones(1,sigE1)*NaN;
peaks_ss=ones(1,sigE2)*NaN;
peaks_us=ones(1,sigE2)*NaN;
vss=ones(sigE2,4)*NaN;
vuss=ones(sigE2,4)*NaN;

% storing all parameters 
data=[ga0 gas gab gb0 gbs gba Jas Jbs Jab Jba nas nbs nab nba gma gmb];

% For pertubing the network: dual signaling to uni-signaling PFL
pa=1.0;             %  pa=0 to kill S->A signal
pb=1.0;             %  pb=0 to kill S->B signal

% Parameter searching begins for bistbale switches
for kk=1:nr         % loop for parameters
    
    % *********************************************************************
    % Finding no. of stable and unstable steady states with variation of
    % bifurcation parameter. Stable and unstable steady states are obtained
    % from the minima and maxima in the potential energy function (PEF)
    % *********************************************************************
    
    % 1-p bifurcation run: low resolution run over singal S
        %******************************************************************
        % if the low-res runs results bistability then a high-res run is 
        % performed and the type of switch is also determined. This allows
        % minimize the computational cost
        %******************************************************************
    for jj=1:sigE1            % loop for bifurcation parameter
        
        S0=sig1(jj)*dsig1;
        AS1=(S0./Jas(kk)).^nas(kk);
        AS_int=1-AS1./(1+AS1);
        BS1=(S0./Jbs(kk)).^nbs(kk);
        BS_int=1-BS1./(1+BS1);
        AB1=(b1./Jab(kk)).^nab(kk);
        AB_int=AB1./(1+AB1);
        Ass=(ga0(kk)+pa.*gas(kk).*AS_int+gab(kk).*AB_int)./gma(kk);
        BA1=(Ass./Jba(kk)).^nba(kk);
        BA_int=BA1./(1+BA1);
        
        % force term in the euqation of B
        f=gb0(kk)+pb.*gbs(kk).*BS_int+gba(kk).*BA_int-gmb(kk).*b1;
        
        % z: potential function
        z=dx1*cumtrapz(-f);
        
        % rescaling of potential w.r.t. the minima of the pontential 
        z1=z-min(z);        % used to determine the unstable steady states         
        z2=-z1;             % used to determine the stable steady states
        
        % storing the number of SS at each signal values
        nmbr_peaks_ss(jj)=numel(findpeaks(z2));     % no. of minima in PEF
        nmbr_peaks_us(jj)=numel(findpeaks(z1));     % no. of maxima in PEF
        
    end
    
    srt=nmbr_peaks_ss(cumsum(nmbr_peaks_ss,2) > 0);
    flipdiff=flip(srt);
    end_zero=flipdiff(cumsum(flipdiff,2) > 0);
    final_diff=flip(end_zero);
    diff_ss=diff(final_diff);
    sum_diffss=sum(abs(diff_ss));

    % *********************************************************************
    % End of low res steady state calculation
    % *********************************************************************
    
    % *********************************************************************
    % Calculation of bifurcation diagram for the parameter combination that
    % generates more than one stable steady states or multistability
    % *********************************************************************
    
    % Irrev. switches and tristability are discarded
    if (nmbr_peaks_ss(1)>1)         % irreversible switch
        continue
    elseif (nmbr_peaks_us(end)>=1)  % irreversible switch
        continue
    elseif (nmbr_peaks_us(1)>=1)    % irreversible switch
        continue
    elseif (nmbr_peaks_ss(end)>1)   % irreversible switch
        continue                    
    elseif any(end_zero==0)
        continue
    elseif max(nmbr_peaks_ss)==3    % tristability    
        disp('TRISTABILITY')
    elseif max(nmbr_peaks_ss)==2    % only bistability is picked up
        kk                         

        % 1-p bifurcation run: high resolution run over singal S
        for j=1:sigE2      % Loop for bifurcation parameter
            
            S0=sig2(j)*dsig2;
            AS1=(S0./Jas(kk)).^nas(kk);
            AS_int=1-AS1./(1+AS1);
            BS1=(S0./Jbs(kk)).^nbs(kk);
            BS_int=1-BS1./(1+BS1);
            AB1=(b2./Jab(kk)).^nab(kk);
            AB_int=AB1./(1+AB1);
            Ass=(ga0(kk)+pa.*gas(kk).*AS_int+gab(kk).*AB_int)./gma(kk);
            BA1=(Ass./Jba(kk)).^nba(kk);
            BA_int=BA1./(1+BA1);
        
            f2=gb0(kk)+pb.*gbs(kk).*BS_int+gba(kk).*BA_int-gmb(kk).*b2;
            
            z0=dx2*cumtrapz(-f2);
            z11=z0-min(z0);
            z22=-z11;
            
            peaks_ss(j)=numel(findpeaks(z22));      % no. of minima in PEF
            peaks_us(j)=numel(findpeaks(z11));      % no. of maxima in PEF
            
            % pks_stb: no. of stable ss
            [pks_stb locs_stb] = findpeaks(z22); 
            
            % storing of stable branches
            if length(locs_stb)==3
                vss(j,1)=S0;
                vss(j,2)=b2(locs_stb(1));
                vss(j,3)=b2(locs_stb(2));
                vss(j,4)=b2(locs_stb(3));
            elseif length(locs_stb)==2
                vss(j,1)=S0;
                vss(j,2)=b2(locs_stb(1));
                vss(j,3)=b2(locs_stb(2));
                vss(j,4)=NaN;
            elseif length(locs_stb)==1
                vss(j,1)=S0;
                vss(j,2)=b2(locs_stb(1));
                vss(j,3)=NaN;
                vss(j,4)=NaN;
            end
            
            % pks_us: no. of unstable ss
            [pks_us locs_us] = findpeaks(z11);   
            
            % storing of unstable branches
            if length(locs_us)==2
                vuss(j,1)=S0;
                vuss(j,2)=b2(locs_us(1));
                vuss(j,3)=b2(locs_us(2));
                vuss(j,4)=NaN;
            elseif length(locs_us)==1
                vuss(j,1)=S0;
                vuss(j,2)=b2(locs_us(1));
                vuss(j,3)=NaN;
                vuss(j,4)=NaN;
            elseif length(locs_us)==0
                vuss(j,1)=S0;
                vuss(j,2)=NaN;
                vuss(j,3)=NaN;
                vuss(j,4)=NaN;
            end
        end
        
        % Discarded: if the high res run genrates tristability 
        if max(peaks_ss)>=3         
            continue
        end

        % saving parameters for all types of multistability in the high-res
        % run
        ic2=ic2+1;
        para_BS_all(ic2,:)=data(kk,:);
        
        % plot of bifurcation diagrams
%         fig1=figure(1);
%         subplot(6,6,ic2)
%         hold off
%         plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%         hold on;
%         plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%         plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
        
        
        % indirect determination of number of regions od bistability
        % sum_diffss2=2 for each bistable region 
        % switches with single BS region: sum_diffss2=2
        % switches with double BS region: sum_diffss2=4
        % switches with tripple BS region: sum_diffss2=6
        srt2=peaks_ss(cumsum(peaks_ss,2) > 0);
        flipd2=flip(srt2);
        rm_zero1=flipd2(cumsum(flipd2,2) > 0);
        rm_zero2=flip(rm_zero1);
        diff_ss2=diff(rm_zero2);
        sum_diffss2=sum(abs(diff_ss2));
        
        % Discarded: if the left most SN point appears before the upper
        % stable branch due to the fixed maximum value of B in the calculation,
        % then the jump in the stable branch at the left most SN point will
        % be missed. Similar situation may appear for the right most SN
        % point. Such cases are discarded
        unF=find(~isnan(vuss(:,2)), 1, 'first');    % first pnt in unS branch
        unL=find(~isnan(vuss(:,2)), 1, 'last');     % last pnt in unS branch
        sbF=find(~isnan(vss(:,2)), 1, 'first');     % first pnt in sS branch
        sbL=find(~isnan(vss(:,2)), 1, 'last');      % last pnt in sS branch
    
        if unF == sbF | unL == sbL 
            continue
        end
        
        % calculation of jump patterens of the stable SS
        vss1=vss;
        vus1=vuss;
        xpeak1=diff(vss1(:,2));
%         figure(6)
%         plot(abs(xpeak1));
%         hold on
        % The fluctuations in the difference (xpeak1) are removed in the BS
        % regions
        i1=0;
        for i1=1:length(xpeak1)
            i2=i1+1;
            if vus1(i2,2)>0 && vus1(i2-1,2)>0
                xpeak1(i1)=0;  
            end
        end
%         figure(6)
%         plot(abs(xpeak1));
         
        xpeak1(abs(xpeak1)<=3)=0;   % small jumps <=3 are removed
        xpeak=xpeak1;
%         figure(6)
%         plot(abs(xpeak));
                
        [Grd2mx Loc_grd2mx]=findpeaks(xpeak);
        [Grd2mn Loc_grd2mn]=findpeaks(-xpeak);
        Grd2mx1=Grd2mx(Grd2mx~=0);
        Grd2mn1=Grd2mn(Grd2mn~=0);
        Loc_grd2mx1=Loc_grd2mx(Grd2mx~=0);
        Loc_grd2mn1=Loc_grd2mn(Grd2mn~=0);
        
        i2=0;
        mxDf=[];
        locmx=[];
        for i1=1:length(Loc_grd2mx1)
            if vus1(Loc_grd2mx1(i1),2)>0
                i2=i2+1;
                mxDf(i2)=Grd2mx1(i1);
                locmx(i2)=Loc_grd2mx1(i1);
            end
        end
        
        i2=0;
        mnDf=[];
        locmn=[];
        for i1=1:length(Loc_grd2mn1)
            if vus1(Loc_grd2mn1(i1)+1,2)>0
                i2=i2+1;
                mnDf(i2)=Grd2mn1(i1);
                locmn(i2)=Loc_grd2mn1(i1);
            end
        end
        
        locmx;
        locmn;
        LoC=[locmx locmn];  % locations of upward and downward jumps
        MxMn=[mxDf mnDf];
        jmpN=length(MxMn);  % number of jupms
        

        % final storing of parameters for bistable switches in the high-res run
        ic1=ic1+1;
        bistability_par(ic1,1)=kk;      % stores the index of selected parameters
        para_BS_pck(ic1,:)=data(kk,:);  % stores the selected parameters
        
        % Plotting the bifurcations segregated according to their type
        

        % seggregation of switches based on the number of bistable regions
        % and jump patterns of stable branch at the SN points
        if (max(peaks_ss)==2 && (sum_diffss2==2))
            if jmpN==0
                disp('Normal Isola')
                iso=iso+1;                      % count of switch
                isopar(iso,1)=kk;               % index for the parameter
                para_iso(iso,:)=data(kk,:);     % parameter for the switch
%                 fig2=figure(2);
%                 subplot(3,4,1);
%                 title('Normal Isola');
%                 plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                 hold on;
%                 plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                 plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
            elseif jmpN==1
                disp('Simple Bistable')
                bis=bis+1;
                bispar(bis,1)=kk;
                para_bis(bis,:)=data(kk,:);
%                 subplot(3,4,2);
%                 title('Simple Bistable');
%                 plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                 hold on;
%                 plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                 plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
            elseif jmpN==2
                disp('Inverted Isola')
                inviso=inviso+1;
                invisopar(inviso,1)=kk;
                para_inviso(inviso,:)=data(kk,:);
%                 subplot(3,4,3);
%                 title('Inverted Isola')
%                 plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                 hold on;
%                 plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                 plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
            end
        elseif (max(peaks_ss)==2 && (sum_diffss2==4))
            if jmpN==0
                disp('Dual Normal Isola')
                dni=dni+1;
                dnipar(dni,1)=kk;
                para_dni(dni,:)=data(kk,:);
%                 subplot(3,4,4);
%                 title('Dual Normal Isola');
%                 plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                 hold on;
%                 plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                 plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
            elseif jmpN==1
                disp('Normal Isola with Bistable')
                nibi=nibi+1;
                nibipar(nibi,1)=kk;
                para_nibi(nibi,:)=data(kk,:);
%                 subplot(3,4,5);
%                 title('NI with Bistable');
%                 plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                 hold on;
%                 plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                 plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
            elseif jmpN==2
                %                 if (loc_grd2mx1(1)<loc_grd2mn1(1) && loc_grd2mx1(2)<loc_grd2mn1(2))
                if length(locmn)==2 | length(locmx)==2
                    disp('Dual Bistable')
                    dbs=dbs+1;
                    dbspar(dbs,1)=kk;
                    para_dbs(dbs,:)=data(kk,:);
%                     subplot(3,4,6);
%                     title('Dual Bistable');
%                     plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                     hold on;
%                     plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                     plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
                elseif locmx(1)>locmn(1)
                    disp('Inverted Mushroom')
                    invmush=invmush+1;
                    invmushpar(invmush,1)=kk;
                    para_invmush(invmush,:)=data(kk,:);
%                     subplot(3,4,7)
%                     title('Inverted Mushroom');
%                     plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                     hold on;
%                     plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                     plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
                elseif locmx(1)<locmn(1)
                    disp('Normal Mushroom')
                    mush=mush+1;
                    mushpar(mush,1)=kk;
                    para_mush(mush,:)=data(kk,:);
%                     subplot(3,4,8);
%                     title('Normal Mushroom');
%                     plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                     hold on;
%                     plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                     plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
                end
            elseif jmpN==3
                disp('Inverted Isola with Bistable')
                iibi=iibi+1;
                iibipar(iibi,1)=kk;
                para_iibi(iibi,:)=data(kk,:);
%                 subplot(3,4,9);
%                 title('II with Bistable');
%                 plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                 hold on;
%                 plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                 plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
            end
        elseif (max(peaks_ss)==2 && (sum_diffss2==6))
            if jmpN==2
                if length(locmx)==1 | length(locmn)==1
                    disp('Normal Isola-Mushroom')
                    nim=nim+1;
                    nimpar(nim,1)=kk;
                    para_nim(nim,:)=data(kk,:);
%                     subplot(3,4,10);
%                     title('NI with Mushroom');
%                     plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                     hold on;
%                     plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                     plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
                elseif length(locmx)==2 | length(locmn)==2
                    disp('Normal Isola-Dual Bistable')
                    nidbs=nidbs+1;
                    nidbspar(nidbs,1)=kk;
                    para_nidbs(nidbs,:)=data(kk,:);
                end     
            elseif jmpN==3
                if length(locmx)==3 | length(locmn)==3
                    disp('Triple Bistable')
                    tbs=tbs+1;
                    tbspar(tbs,1)=kk;
                    para_tbs(tbs,:)=data(kk,:);
                else
                    disp('Mushroom with Bistable')
                    mbs=mbs+1;
                    mbspar(mbs,1)=kk;
                    para_mbs(mbs,:)=data(kk,:);
%                     subplot(3,4,11);
%                     title('Mushroom with Bistable');
%                     plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                     hold on;
%                     plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                     plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
                end
            elseif jmpN==4
                if length(locmx)==2 | length(locmn)==2
                    disp('Inverted Isola-Mushroom')
                    iim=iim+1;
                    iimpar(iim,1)=kk;
                    para_iim(iim,:)=data(kk,:);
%                     subplot(3,4,12);
%                     title('II with Mushroom');
%                     plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
%                     hold on;
%                     plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
%                     plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
                elseif length(locmx)==3 | length(locmn)==3
                   disp('Inverted Isola-Dual Bistable')
                   iidbs=iidbs+1;
                   iidbspar(iidbs,1)=kk;
                   para_iidbs(iidbs,:)=data(kk,:);
                end
            end
        end
        
        clear stbss unsss vss vuss vss1 vssN vus1 Grd2mn Grd2mn1 Grd2mx ...
            Grd2mx1 Loc_grd2mn Loc_grd2mn1 Loc_grd2mx Loc_grd2mx1 ...
            mxDf mnDf locmn locmx LoC MxMn
        
    end
end
% set(fig1,'Visible','on');
% set(fig2,'Visible','on');

tot_bis=bis;
tot_iso=iso;
tot_inviso=inviso;
tot_dni=dni;
tot_nibi=nibi;
tot_dbs=dbs;
tot_mush=mush;
tot_invmush=invmush;
tot_iibi=iibi;
tot_nim=nim;
tot_mbs=mbs;
tot_iim=iim;
tot_nidbs=nidbs;
tot_tbs=tbs;
tot_iidbs=iidbs;

% seggregated counts of different types of switches
nnMA_OR_dataR.total_bistability=[tot_bis tot_iso tot_inviso tot_dni tot_nibi ...
    tot_dbs tot_mush tot_invmush tot_iibi tot_nim tot_mbs tot_iim tot_nidbs ...
      tot_tbs tot_iidbs];           
total_seg=sum(nnMA_OR_dataR.total_bistability)
% no. of discarded switches due to tristability and assymptotic rise of B
discrepency=total_seg-length(bistability_par) 

% storing of seggregated parameters
nnMA_OR_dataR.data=data;
nnMA_OR_dataR.para_BS_all=para_BS_all;
nnMA_OR_dataR.para_BS_pck=para_BS_pck;
nnMA_OR_dataR.para_bis=para_bis;
nnMA_OR_dataR.para_iso=para_iso;
nnMA_OR_dataR.para_invis=para_inviso;
nnMA_OR_dataR.para_dni=para_dni;
nnMA_OR_dataR.para_nibi=para_nibi;
nnMA_OR_dataR.para_dbs=para_dbs;
nnMA_OR_dataR.para_mush=para_mush;
nnMA_OR_dataR.para_invmush=para_invmush;
nnMA_OR_dataR.para_iibi=para_iibi;
nnMA_OR_dataR.para_nim=para_nim;
nnMA_OR_dataR.para_mbs=para_mbs;
nnMA_OR_dataR.para_iim=para_iim;
nnMA_OR_dataR.para_nidbs=para_nidbs;
nnMA_OR_dataR.para_tbs=para_tbs;
nnMA_OR_dataR.para_iidbs=para_iidbs;

% storing of index of seggregated parameters
nnMA_OR_dataR.indx_bis=bispar;
nnMA_OR_dataR.indx_iso=isopar;
nnMA_OR_dataR.indx_invis=invisopar;
nnMA_OR_dataR.indx_dni=dnipar;
nnMA_OR_dataR.indx_nibi=nibipar;
nnMA_OR_dataR.indx_dbs=dbspar;
nnMA_OR_dataR.indx_mush=mushpar;
nnMA_OR_dataR.indx_invmush=invmushpar;
nnMA_OR_dataR.indx_iibi=iibipar;
nnMA_OR_dataR.indx_nim=nimpar;
nnMA_OR_dataR.indx_mbs=mbspar;
nnMA_OR_dataR.indx_iim=iimpar;
nnMA_OR_dataR.indx_nidbs=nidbspar;
nnMA_OR_dataR.indx_tbs=tbspar;
nnMA_OR_dataR.indx_iidbs=iidbspar;

save parscan_nnMA_OR_R1.mat nnMA_OR_dataR

toc

% quit
