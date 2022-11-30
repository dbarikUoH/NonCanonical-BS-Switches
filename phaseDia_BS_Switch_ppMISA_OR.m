% Bifurcation analysis of feedback refulated systems using potential
% energy function (PEF) calculation.

% Phase diagrams of canonical and noncanonical bistable switches

% Date: November 25, 2022
% Debashis Barik, School of Chemistry, University of Hyderabad, India

% Model: ppMISA_OR
% Logic Gate: OR

clear
clc
tic

% Initial values of the counts of reversible and irreversible switches
iso=0;              % Switch: Is                            Sr. No.: 1
bsF=0;              % Switch: BsF                           Sr. No.: 2
bsB=0;              % Switch: BsB                           Sr. No.: 3
inis=0;             % Switch: InIs                          Sr. No.: 4
dlis=0;             % Switch: DIs                           Sr. No.: 5
dlinis=0;           % Switch: DInIs                         Sr. No.: 6
dlbsF=0;            % Switch: DBsF                          Sr. No.: 7
dlbsB=0;            % Switch: DBsB                          Sr. No.: 8
mus=0;              % Switch: Msh                           Sr. No.: 9
inmus=0;            % Switch: InMsh                         Sr. No.: 10    
isbsFb=0;           % Switch: Is-BsF                        Sr. No.: 11
isbsFa=0;           % Switch: BsF-Is                        Sr. No.: 11
inisbsFb=0;         % Switch: InIs-BsF                      Sr. No.: 12
inisbsFa=0;         % Switch: BsF-InIs                      Sr. No.: 13
isbsBb=0;           % Switch: Is-BsB                        Sr. No.: 14
inisbsBb=0;         % Switch: InIs-BsB                      Sr. No.: 15
isbsBa=0;           % Switch: BsB-Is                        Sr. No.: 16
inisbsBa=0;         % Switch: BsB-InIs                      Sr. No.: 17
tpbsF=0;            % Switch: TBsF                          Sr. No.: 18
tpbsB=0;            % Switch: TBsB                          Sr. No.: 19
isdlbsFb=0;         % Switch: Is-DBsF                       Sr. No.: 20
isdlbsFm=0;         % Switch: BsF-Is-BsF                    Sr. No.: 21
isdlbsFa=0;         % Switch: DBsF-Is                       Sr. No.: 22    
inisdlbsFb=0;       % Switch: InIs-DBsF                     Sr. No.: 23    
inisdlbsFm=0;       % Switch: BsF-InIs-BsF                  Sr. No.: 24
inisdlbsFa=0;       % Switch: DBsF-InIs                     Sr. No.: 25
isdlbsBb=0;         % Switch: Is-DBsB                       Sr. No.: 26
isdlbsBm=0;         % Switch: BsB-Is-BsB                    Sr. No.: 27
isdlbsBa=0;         % Switch: DBsB-Is                       Sr. No.: 28  
inisdlbsBb=0;       % Switch: InIs-DBsB                     Sr. No.: 29    
inisdlbsBm=0;       % Switch: BsB-InIs-BsB                  Sr. No.: 30
inisdlbsBa=0;       % Switch: DBsB-InIs                     Sr. No.: 31    
dlisbsFbb=0;        % Switch: DIs-BsF                       Sr. No.: 32
dlisbsFba=0;        % Switch: Is-BsF-Is                     Sr. No.: 33
dlisbsFaa=0;        % Switch: BsF-DIs                       Sr. No.: 34
dlisbsBbb=0;        % Switch: DIs-BsB                       Sr. No.: 35
dlisbsBba=0;        % Switch: Is-BsB-Is                     Sr. No.: 36
dlisbsBaa=0;        % Switch: BsB-DIs                       Sr. No.: 37
ismusb=0;           % Switch: Is-Msh                        Sr. No.: 38
ismusm=0;           % Switch: BsF-Is-BsB                    Sr. No.: 39
ismusa=0;           % Switch: Msh-Is                        Sr. No.: 40
inismusb=0;         % Switch: InIs-Msh                      Sr. No.: 41
inismusm=0;         % Switch: BsF-InIs-BsB                  Sr. No.: 42
inismusa=0;         % Switch: Msh-InIs                      Sr. No.: 43
isinmusb=0;         % Switch: Is-InMsh                      Sr. No.: 44    
isinmusm=0;         % Switch: BsB-Is-BsF                    Sr. No.: 45
isinmusa=0;         % Switch: InMsh-Is                      Sr. No.: 46
inisinmusb=0;       % Switch: InIs-InMsh                    Sr. No.: 47    
inisinmusm=0;       % Switch: BsB-InIs-BsF                  Sr. No.: 48
inisinmusa=0;       % Switch: InMsh-InIs                    Sr. No.: 49
musbsFa=0;          % Switch: Msh-BsF                       Sr. No.: 50
musbsBa=0;          % Switch: Msh-BsB                       Sr. No.: 51
musbsFb=0;          % Switch: InMsh-BsF                     Sr. No.: 52
inmusbsBa=0;        % Switch: BsF-Msh                       Sr. No.: 53
inmusbsFa=0;        % Switch: InMsh-BsB                     Sr. No.: 54
inmusbsBb=0;        % Switch: BsB-InMsh                     Sr. No.: 55
tlis=0;             % Switch: Tis                           Sr. No.: 56
LbsB=0;             % Switch: L-BsB                         Sr. No.: 57
LbsF=0;             % Switch: L-BsF                         Sr. No.: 58
Ldlis=0;            % Switch: L-DIs                         Sr. No.: 59
LbsFinis=0;         % Switch: L-DInIs                       Sr. No.: 60
LdlbsF=0;           % Switch: L-DBsF                        Sr. No.: 61
LbsBis=0;           % Switch: L-DBsB                        Sr. No.: 62
Lmus=0;             % Switch: L-Msh                         Sr. No.: 63
Linmus=0;           % Switch: L-InMsh                       Sr. No.: 64
LbsFis=0;           % Switch: L-BsF-Is                      Sr. No.: 65
LbsBinis=0;         % Switch: L-BsB-InIs                    Sr. No.: 66
LtpbsB=0;           % Switch: L-TBsB                        Sr. No.: 67
LinmusbsFa=0;       % Switch: L-Is-DBsF                     Sr. No.: 68
LdlbsFism=0;        % Switch: L-BsF-Is-BsF                  Sr. No.: 69
LdlbsFisa=0;        % Switch: L-DBsF-Is                     Sr. No.: 70
LdlbsBism=0;        % Switch: L-BsB-Is-BsB                  Sr. No.: 71
LdlbsBisa=0;        % Switch: L-DBsB-Is                     Sr. No.: 72
Linmusism=0;        % Switch: L-Dis-BsF                     Sr. No.: 73
Linmusisa=0;        % Switch: L-Is-BsF-Is                   Sr. No.: 74
LdlisbsF=0;         % Switch: L-BsF-DIs                     Sr. No.: 75
LdlisbsB=0;         % Switch: L-BsB-DIs                     Sr. No.: 76
LmusbsB=0;          % Switch: L-Is-Msh                      Sr. No.: 77
Lmusism=0;          % Switch: L-BsF-Is-BsB                  Sr. No.: 78
Lmusisa=0;          % Switch: L-Msh-Is                      Sr. No.: 79
LbsBinmus=0;        % Switch: L-Is-InMsh                    Sr. No.: 80
RbsF=0;             % Switch: BsF-R                         Sr. No.: 81
RbsB=0;             % Switch: BsB-R                         Sr. No.: 82
Rdlis=0;            % Switch: R-DIs                         Sr. No.: 83
Rdlinis=0;          % Switch: DInIs-R                       Sr. No.: 84
RbsFis=0;           % Switch: DBsF-R                        Sr. No.: 85
RbsBinis=0;         % Switch: DBsB-R                        Sr. No.: 86
Rmus=0;             % Switch: Msh-R                         Sr. No.: 87
Rinmus=0;           % Switch: InMsh-R                       Sr. No.: 88
RbsFinis=0;         % Switch: InIs-BsF-R                    Sr. No.: 89
RbsBis=0;           % Switch: Is-BsB-R                      Sr. No.: 90
RtpbsF=0;           % Switch: TBsF-R                        Sr. No.: 91  
RdlbsFisb=0;        % Switch: Is-DBsF-R                     Sr. No.: 92
RdlbsFism=0;        % Switch: BsF-Is-BsF-R                  Sr. No.: 93
RdlbsBisb=0;        % Switch: Is-DBsB-R                     Sr. No.: 94
RdlbsBism=0;        % Switch: BsB-Is-BsB-R                  Sr. No.: 95  
RinmusbsB=0;        % Switch: DBsB-Is-R                     Sr. No.: 96    
RdlisbsF=0;         % Switch: DIs-BsF-R                     Sr. No.: 97
RdlisbsB=0;         % Switch: Dis-BsB-R                     Sr. No.: 98
Rinmusisb=0;        % Switch: Is-BsB-Is-R                   Sr. No.: 99
Rinmusism=0;        % Switch: BsB-DIs-R                     Sr. No.: 100
Rmusisb=0;          % Switch: Is-Msh-R                      Sr. No.: 101
Rmusism=0;          % Switch: BsF-Is-BsB-R                  Sr. No.: 102
RmusbsF=0;          % Switch: Msh-Is-R                      Sr. No.: 103
RinmusbsF=0;        % Switch: InMsh-Is-R                    Sr. No.: 104
LRbs=0;             % Switch: L-Bs-R                        Sr. No.: 105
LRisis=0;           % Switch: L-DIs-R                       Sr. No.: 106
LRinisinis=0;       % Switch: L-DInIs-R                     Sr. No.: 107
LRinisis=0;         % Switch: L-DBsF-R                      Sr. No.: 108
LRisinis=0;         % Switch: L-DBsB-R                      Sr. No.: 109
LRinisbsFis=0;      % Switch: L-TBsF-R                      Sr. No.: 110
LRisbsBinis=0;      % Switch: L-TBsB-R                      Sr. No.: 111    
LRisbsFis=0;        % Switch: L-Is-DBsF-R                   Sr. No.: 112
LRinisisis=0;       % Switch: L-BsF-Is-BsF-R                Sr. No.: 113
LRinisinisis=0;     % Switch: L-BsF-InIs-BsF-R              Sr. No.: 114
LRinisbsFinis=0;    % Switch: L-DBsF-InIs-R                 Sr. No.: 115
LRisisinis=0;       % Switch: L-BsB-Is-BsB-R                Sr. No.: 116
LRisbsBis=0;        % Switch: L-DBsB-Is-R                   Sr. No.: 117
LRinisbsBinis=0;    % Switch: L-InIs-DBsB-R                 Sr. No.: 118    
LRisinisinis=0;     % Switch: L-BsB-InIs-BsB-R              Sr. No.: 119
LRisisis=0;         % Switch: L-DIs-BsF-R                   Sr. No.: 120
LRisbsFinis=0;      % Switch: L-Is-Msh-R                    Sr. No.: 122
LRinisisinis=0;     % Switch: L-BsF-Is-BsB-R                Sr. No.: 123
LRinisbsBis=0;      % Switch: L-Msh-Is-R                    Sr. No.: 124
LRinisinisinis=0;   % Switch: L-BsF-InIs-BsB-R              Sr. No.: 125
LRisinisis=0;       % Switch: L-BsB-InIs-BsF-R              Sr. No.: 126

% These four ar extra: was not counted in the paper
LbsBisinism=0;
LbsBisinisa=0;
RisinisbsF=0;
RinisisbsF=0;

% Array for searched parameters for differen switches
para_BS_all=[];
para_BS_pck=[];
para_iso=[];
para_bsF=[];
para_bsB=[];
para_inis=[];
para_dlis=[];
para_isbsBa=[];
para_isbsFa=[];
para_isbsBb=[];
para_isbsFb=[];
para_dlbsF=[];
para_dlbsB=[];
para_mus=[];
para_inmus=[];
para_inisbsFb=[];
para_inisbsFa=[];
para_inisbsBb=[];
para_inisbsBa=[];
para_dlinis=[];
para_tlis=[];
para_dlisbsFbb=[];
para_dlisbsFba=[];
para_dlisbsFaa=[];
para_dlisbsBbb=[];
para_dlisbsBba=[];
para_dlisbsBaa=[];
para_ismusb=[];
para_ismusm=[];
para_ismusa=[];
para_isdlbsFb=[];
para_isdlbsFm=[];
para_isdlbsFa=[];
para_isdlbsBb=[];
para_isdlbsBm=[];
para_isdlbsBa=[];
para_isinmusb=[];
para_isinmusm=[];
para_isinmusa=[];
para_tpbsF=[];
para_tpbsB=[];
para_musbsFa=[];
para_musbsBa=[];
para_musbsFb=[];
para_inmusbsBa=[];
para_inmusbsFa=[];
para_inmusbsBb=[];
para_inisinmusb=[];
para_inisinmusm=[];
para_inisinmusa=[];
para_inismusb=[];
para_inismusm=[];
para_inismusa=[];
para_inisdlbsFb=[];
para_inisdlbsFm=[];
para_inisdlbsFa=[];
para_inisdlbsBb=[];
para_inisdlbsBm=[];
para_inisdlbsBa=[];
para_LbsB=[];
para_LbsF=[];
para_Ldlis=[];
para_LbsFis=[];
para_LbsBis=[];
para_Linmus=[];
para_LdlbsF=[];
para_Lmus=[];
para_LbsBinis=[];
para_LbsFinis=[];
para_RbsF=[];
para_RbsB=[];
para_Rdlis=[];
para_RbsBis=[];
para_RbsFis=[];
para_Rinmus=[];
para_RbsFinis=[];
para_Rmus=[];
para_RbsBinis=[];
para_Rdlinis=[];
para_LdlisbsB=[];
para_LdlisbsF=[];
para_LdlbsBisa=[];
para_Linmusisa=[];
para_LdlbsBism=[];
para_Linmusism=[];
para_Lmusisa=[];
para_LdlbsFisa=[];
para_Lmusism=[];
para_LdlbsFism=[];
para_LbsBisinism=[];
para_LtpbsB=[];
para_LbsBinmus=[];
para_LmusbsB=[];
para_LinmusbsFa=[];
para_LbsBisinisa=[];
para_RdlisbsF=[];
para_Rinmusism=[];
para_RdlbsFism=[];
para_Rinmusisb=[];
para_RdlbsFisb=[];
para_RdlisbsB=[];
para_RinisisbsF=[];
para_RinmusbsB=[];
para_RinmusbsF=[];
para_RdlbsBism=[];
para_RmusbsF=[];
para_RtpbsF=[];
para_Rmusism=[];
para_RisinisbsF=[];
para_RdlbsBisb=[];
para_Rmusisb=[];
para_LRbs=[];
para_LRisis=[];
para_LRinisis=[];
para_LRisinis=[];
para_LRinisinis=[];
para_LRisisis=[];
para_LRinisisis=[];
para_LRisbsBis=[];
para_LRisbsFis=[];
para_LRisisinis=[];
para_LRinisbsBis=[];
para_LRinisbsFis=[];
para_LRinisisinis=[];
para_LRisinisis=[];
para_LRisbsBinis=[];
para_LRisbsFinis=[];
para_LRinisinisis=[];
para_LRinisbsBinis=[];
para_LRinisbsFinis=[];
para_LRisinisinis=[];
para_LRinisinisinis=[];

iso_par=[];
bsF_par=[];
bsB_par=[];
inis_par=[];
dlis_par=[];
isbsBa_par=[];
isbsFa_par=[];
isbsBb_par=[];
isbsFb_par=[];
dlbsF_par=[];
dlbsB_par=[];
mus_par=[];
inmus_par=[];
inisbsFb_par=[];
inisbsFa_par=[];
inisbsBb_par=[];
inisbsBa_par=[];
dlinis_par=[];
tlis_par=[];
dlisbsFbb_par=[];
dlisbsFba_par=[];
dlisbsFaa_par=[];
dlisbsBbb_par=[];
dlisbsBba_par=[];
dlisbsBaa_par=[];
ismusb_par=[];
ismusm_par=[];
ismusa_par=[];
isdlbsFb_par=[];
isdlbsFm_par=[];
isdlbsFa_par=[];
isdlbsBb_par=[];
isdlbsBm_par=[];
isdlbsBa_par=[];
isinmusb_par=[];
isinmusm_par=[];
isinmusa_par=[];
tpbsF_par=[];
tpbsB_par=[];
musbsFa_par=[];
musbsBa_par=[];
musbsFb_par=[];
inmusbsBa_par=[];
inmusbsFa_par=[];
inmusbsBb_par=[];
inisinmusb_par=[];
inisinmusm_par=[];
inisinmusa_par=[];
inismusb_par=[];
inismusm_par=[];
inisdlbsFb_par=[];
inisdlbsFm_par=[];
inisdlbsFa_par=[];
inisdlbsBb_par=[];
inisdlbsBm_par=[];
inisdlbsBa_par=[];
LbsB_par=[];
LbsF_par=[];
Ldlis_par=[];
LbsFis_par=[];
LbsBis_par=[];
Linmus_par=[];
LdlbsF_par=[];
Lmus_par=[];
LbsBinis_par=[];
LbsFinis_par=[];
RbsF_par=[];
RbsB_par=[];
Rdlis_par=[];
RbsBis_par=[];
RbsFis_par=[];
Rinmus_par=[];
RbsFinis_par=[];
Rmus_par=[];
RbsBinis_par=[];
Rdlinis_par=[];
LdlisbsB_par=[];
LdlisbsF_par=[];
LdlbsBisa_par=[];
Linmusisa_par=[];
LdlbsBism_par=[];
Linmusism_par=[];
Lmusisa_par=[];
LdlbsFisa_par=[];
Lmusism_par=[];
LdlbsFism_par=[];
LbsBisinism_par=[];
LtpbsB_par=[];
LbsBinmus_par=[];
LmusbsB_par=[];
LinmusbsFa_par=[];
LbsBisinisa_par=[];
RdlisbsF_par=[];
Rinmusism_par=[];
RdlbsFism_par=[];
Rinmusisb_par=[];
RdlbsFisb_par=[];
RdlisbsB_par=[];
RinisisbsF_par=[];
RinmusbsB_par=[];
RinmusbsF_par=[];
RdlbsBism_par=[];
RmusbsF_par=[];
RtpbsF_par=[];
Rmusism_par=[];
RisinisbsF_par=[];
RdlbsBisb_par=[];
Rmusisb_par=[];
LRbs_par=[];
LRisis_par=[];
LRinisis_par=[];
LRisinis_par=[];
LRinisinis_par=[];
LRisisis_par=[];
LRinisisis_par=[];
LRisbsBis_par=[];
LRisbsFis_par=[];
LRisisinis_par=[];
LRinisbsBis_par=[];
LRinisbsFis_par=[];
LRinisisinis_par=[];
LRisinisis_par=[];
LRisbsBinis_par=[];
LRisbsFinis_par=[];
LRinisinisis_par=[];
LRinisbsBinis_par=[];
LRinisbsFinis_par=[];
LRisinisinis_par=[];
LRinisinisinis_par=[];

bistability_par=[];

ic1=0;
ic2=0;

rng('shuffle')

load parscan_ppMISA_OR_R1.mat
parM=ppMISA_OR_dataR.para_mush;

% parM=ppMI_OR_dataR.para_mush and kk=7 for phase dia of isola
% parM=ppMI_OR_dataR.para_invmush and kk=12 for phase dia of inv-isola

nr=size(parM,1);                    % no. of parameter combinations

sigE1=4000;
sigE2=20000;
dsig1=0.25;
dsig2=0.05;
dx1=1;
dx2=1;
xE=10000;
b1=[0:dx1:xE];              % range for argument in potential function
b2=[0:dx2:xE];

ss=[];
uss=[];
sig1=[1:sigE1];
sig2=[1:sigE2];

gba0=15;
gab0=30;
dp=5;
dim1=20;
dim2=20;
phaseIndx=ones(dim1,dim2)*NaN;
Val_gba=ones(dim1,dim2)*NaN;
Val_gab=ones(dim1,dim2)*NaN;

for kk=1:1               % loop for parameters
    kk=7
    ga0=parM(kk,1);
    gas=parM(kk,2);
    gab=parM(kk,3);
    gb0=parM(kk,4);
    gbs=parM(kk,5);
    gba=parM(kk,6);
    gbb=parM(kk,7);
    Jas=parM(kk,8);
    Jbs=parM(kk,9);
    Jab=parM(kk,10);
    Jba=parM(kk,11);
    Jbb=parM(kk,12);
    nas=parM(kk,13);
    nbs=parM(kk,14);
    nab=parM(kk,15);
    nba=parM(kk,16);
    nbb=parM(kk,17);
    gma=parM(kk,18);
    gmb=parM(kk,19);
    
    fcnt=0;
    for kp1=1:dim1
        gba=gba0+(kp1-1)*dp;
        fnt1=dim1*(dim2-kp1);
        for kp2=1:dim2
            [kp1 kp2]
            fnt2=fnt1+kp2;
            %**************************************************************
            % for manual choice of gab, gba
            %**************************************************************
%              gab=32;
%              gba=90;
            %**************************************************************
            gab=gab0+(kp2-1)*dp;
            Val_gba(kp1,kp2)=gba;
            Val_gab(kp1,kp2)=gab;
            fcnt=fcnt+1;
            nmbr_peaks_ss=ones(1,sigE1)*NaN;
            nmbr_peaks_us=ones(1,sigE1)*NaN;
            dsig2=0.05;
            
            PotEr=[];
            xss=ones(length(sigE1),3)*NaN;
            xus=ones(length(sigE1),3)*NaN;
            % data=[ga0 gas gab gb0 gbs gba gbb Jas Jbs Jab Jba Jbb nas nbs nab nba nbb gma gmb];
            
            % *********************************************************************
            % Finding no. of stable OR unstable steady states with variation of
            % bifurcation parameter. Stable OR unstable steady states are given
            % by minima OR maxima in the PEF respectively.
            % *********************************************************************
            for jj=1:sigE1            % loop for bifurcation parameter
                
                S0=sig1(jj)*dsig1;
                AS1=(S0./Jas).^nas;
                AS_int=AS1./(1+AS1);
                BS1=(S0./Jbs).^nbs;
                BS_int=BS1./(1+BS1);
                AB1=(b1./Jab).^nab;
                AB_int=1-AB1./(1+AB1);
                Ass=(ga0+gas.*AS_int+gab.*AB_int)./gma;
                BA1=(Ass./Jba).^nba;
                BA_int=1-BA1./(1+BA1);
                BB1=(b1./Jbb).^nbb;
                BB_int=BB1./(1+BB1);
                
                f = gb0+gbs.*BS_int+gba.*BA_int+gbb.*BB_int-gmb.*b1;
                
                z=dx1*cumtrapz(-f);   % z: potential function OR it is
                % the negative of the function (f) integrated
                z1=z-min(z);          % rescaling of potential against the global
                % minima of the pontential function
                z2=-z1;
                % figure(2)
                % plot(x,z2)
                % findpeaks(z2)
                nmbr_peaks_ss(jj)=numel(findpeaks(z2));     % no. of minima in PEF
                nmbr_peaks_us(jj)=numel(findpeaks(z1));     % no. of mamima in PEF
                [stb_pks stb_locs] = findpeaks(z2);
                [ustb_pks ustb_locs] = findpeaks(z1);
                if length(stb_locs)==2
                    xss(jj,1)=S0;
                    xss(jj,2)=b1(stb_locs(1));
                    xss(jj,3)=b1(stb_locs(2));
                elseif length(stb_locs)==1
                    xss(jj,1)=S0;
                    xss(jj,2)=b1(stb_locs(1));
                    xss(jj,3)=NaN;
                end
                if length(ustb_locs)==1
                    xus(jj,1)=S0;
                    xus(jj,2)=b1(ustb_locs(1));
                    xus(jj,3)=NaN;
                elseif length(ustb_locs)==0
                    xus(jj,1)=S0;
                    xus(jj,2)=NaN;
                    xus(jj,3)=NaN;
                end
                
            end
            
            srt=nmbr_peaks_ss(cumsum(nmbr_peaks_ss,2) > 0);
            flipdiff=flip(srt);
            end_zero=flipdiff(cumsum(flipdiff,2) > 0);
            final_diff=flip(end_zero);
            diff_ss=diff(final_diff);
            sum_diffss=sum(abs(diff_ss));
            
            %     xx=(nmbr_peaks_ss(~isnan(nmbr_peaks_ss)));
            %     xx1=sum(abs(diff(xx)));
            %     [sum_diffss xx1]
            
            %   Determination of maximum X_ss based on the sparsed calculation
            xmax=max(max(xss(:,2)),max(xss(:,3)));
            xE1=round(xmax,-3)+500;
            b2=[0:dx2:xE1];         % Setting new range of X based on X_ss maximum
            
            %     sum_ssdiff(kk,1)=sum_diffss;
            % *********************************************************************
            % End of finding no. of steady states calculation
            % *********************************************************************
            
            % *********************************************************************
            % Calculation of bifurcation diagram for the parameter combination that
            % generates more than one stable steady states(only for multistability)
            % *********************************************************************
            
            % Counting the no. of stable steady state
            if max(nmbr_peaks_ss)==1
                phaseIndx(kp1,kp2)=0;
            elseif max(nmbr_peaks_ss)==2
                kk;                         % that results bistability
                kk0=find(nmbr_peaks_ss==2,1,'first');
                kk1=find(nmbr_peaks_ss==2,1,'last');
                dS=abs(kk0-kk1);
                if (dS <= 50)
                    sCut=15;
                elseif (50 < dS && dS <= 100)
                    sCut=25;
                elseif (100 < dS && dS <= 200)
                    sCut=35;
                elseif (dS > 200)
                    sCut=45;
                end
                %         kk3=0.9*kk0;
                kk3=kk0-sCut;
                %         kk4=1.1*kk1;
                kk4=kk1+sCut;
                kkF=round(kk3*(sigE2/sigE1),0);
                kkL=round(kk4*(sigE2/sigE1),0);
                sigE3F=max(kkF,1);
                sigE3=min(kkL,sigE2);
                %         sigE3=sigE2;
                
                if sum_diffss == 2
                    if (dS <= 5)            % Delta S=1.25
                        dsig2nw=0.01;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    elseif (5 < dS && dS <=50)
                        dsig2nw=0.025;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    elseif (50 < dS && dS <= 100)
                        dsig2nw=0.05;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    elseif (100 < dS && dS <= 400)
                        dsig2nw=0.1;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    else
                        dsig2nw=0.2;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    end
                elseif sum_diffss ~ 2;
                    if (dS <= 5)
                        dsig2nw=0.01;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    elseif (5<dS && dS<=50)
                        dsig2nw=0.025;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    elseif (50 < dS && dS <= 100)
                        dsig2nw=0.05;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    else
                        dsig2nw=0.1;
                        factr=dsig2nw/dsig2;
                        dsig2=dsig2nw;
                        sigE3=round(sigE3/factr,0);
                        sigE3F1=round(sigE3F/factr,0);
                    end
                elseif sum_diffss == 0;
                    dsig2nw=10;
                    factr=dsig2nw/dsig2;
                    dsig2=dsig2nw;
                    sigE3=round(sigE3/factr,0);
                    sigE3F1=round(sigE3F/factr,0);
                end
                [sum_diffss dS dsig2nw];
                
                if(sigE3F==1)
                    sigE3F1=sigE3F;
                end
                %**********************************************************
                % For manually changing the S region
                %**********************************************************
%                  sigE3F1=800
%                  sigE3=3300
                %**********************************************************
                sZ=sigE3-sigE3F1;
                peaks_ss=ones(1,sZ)*NaN;
                peaks_us=ones(1,sZ)*NaN;
                vss=ones(sZ,4)*NaN;
                vuss=ones(sZ,4)*NaN;
                sig2=[1:sigE3];
                
                
                for j=1:sZ      % Loop for bifurcation parameter
                    
                    j1=sigE3F1+(j-1);
                    S0=sig2(j1)*dsig2;
                    AS1=(S0./Jas).^nas;
                    AS_int=AS1./(1+AS1);
                    BS1=(S0./Jbs).^nbs;
                    BS_int=BS1./(1+BS1);
                    AB1=(b2./Jab).^nab;
                    AB_int=1-AB1./(1+AB1);
                    Ass=(ga0+gas.*AS_int+gab.*AB_int)./gma;
                    BA1=(Ass./Jba).^nba;
                    BA_int=1-BA1./(1+BA1);
                    BB1=(b2./Jbb).^nbb;
                    BB_int=BB1./(1+BB1);
                    
                    f2 = gb0+gbs.*BS_int+gba.*BA_int+gbb.*BB_int-gmb.*b2;
                    
                    z0=dx2*cumtrapz(-f2);
                    z11=z0-min(z0);
                    z22=-z11;
                    
                    %             PotEr(j,:)=z11;                      % potenital energy surface
                    
                    peaks_ss(j)=numel(findpeaks(z22));     % no. of minima in PEF
                    peaks_us(j)=numel(findpeaks(z11));
                    
                    [pks_stb locs_stb] = findpeaks(z22); % pks_stb: no. of stable ss
                    
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
                    
                    [pks_us locs_us] = findpeaks(z11);   % pks_stb: no. of unstable ss
                    
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
                
                if max(peaks_ss)>=3         % to consider only bistable switches
                    continue
                end
                
                fig1=figure(1);
                subplot(dim1,dim2,fnt2)
                %         hold off
                %         contour(log(PotEr'),'LevelStep',0.4,'Linewidth',0.75)
                %         set(gca,'TickLength',[0.03 0.025],'Clim',[-9,9])
                hold on
                plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
                hold on;
                plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
                plot(vss(:,1),vss(:,4),'.','color','k','MarkerSize',4);
                plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
                plot(vuss(:,1),vuss(:,3),'.','color','r','MarkerSize',4);
                set(gca,'XtickLabel',[])
                set(gca,'YtickLabel',[])
                %         set(fig1,'visible','off')
                
                %       saving data for all bistability
                ic2=ic2+1;
                para_BS_all(ic2,:)=parM(kk,:);
                
                srt2=peaks_ss(cumsum(peaks_ss,2) > 0);
                flipd2=flip(srt2);
                rm_zero1=flipd2(cumsum(flipd2,2) > 0);
                rm_zero2=flip(rm_zero1);
                diff_ss2=diff(rm_zero2);
                sum_diffss2=sum(abs(diff_ss2));
                %         xx=(peaks_ss(~isnan(peaks_ss)));
                %         xx1=sum(abs(diff(xx)));
                %         [sum_diffss2 xx1]
                
                % Checking that the stable branch starts before the unstable branch
                % if not then it is discarded from the segragationa and count
                unF=find(~isnan(vuss(:,2)), 1, 'first');
                unL=find(~isnan(vuss(:,2)), 1, 'last');
                sbF=find(~isnan(vss(:,2)), 1, 'first');
                sbL=find(~isnan(vss(:,2)), 1, 'last');
                %        figure(1)
                %        plot(vuss(unF:unL,1),vuss(unF:unL,2))
                if (dsig2nw==0.2)
                    nCut=10;
                elseif (dsig2nw==0.1)
                    nCut=20;
                elseif (dsig2nw==0.05)
                    nCut=30;
                elseif (dsig2nw==0.025)
                    nCut=40;
                else
                    nCut=50;
                end
                nF=max(unF-nCut,1);
                nL=min(unL+nCut,sZ);
                nSize=nL-nF+1;
                %         if unF == sbF | unL == sbL
                %             continue
                %         end
                
                vss1=vss(nF:nL,:);
                vus1=vuss(nF:nL,:);
                
                xpeak1=diff(vss1(:,2));
                %         figure(6)
                %         plot(abs(xpeak1));
                %         hold on
                % The fluctuations in the difference (xpeak2) are removed in the BS regions
                i1=0;
                for i1=1:length(xpeak1)
                    i2=i1+1;
                    if vus1(i2,2)>0 && vus1(i2-1,2)>0
                        xpeak1(i1)=0;
                    end
                end
                %         figure(6)
                %         plot(abs(xpeak1));
                
                xpeak1(abs(xpeak1)<=4)=0;   % small jumps <=4 are removed
                xpeak=xpeak1;
                %         figure(6)
                %         plot(abs(xpeak));
                
                %         [Grd2mx Loc_grd2mx]=findpeaks(abs(xpeak));
                [Grd2mx Loc_grd2mx]=findpeaks(xpeak);
                [Grd2mn Loc_grd2mn]=findpeaks(-xpeak);
                Grd2mx1=Grd2mx(Grd2mx~=0);
                Grd2mn1=Grd2mn(Grd2mn~=0);
                Loc_grd2mx1=Loc_grd2mx(Grd2mx~=0)+1;    % 1 is added to correct the location
                Loc_grd2mn1=Loc_grd2mn(Grd2mn~=0)+1;    % 1 is added to correct the location
                %         Grd2mn1=[];
                %         Loc_grd2mn1=[];
                
                i2=0;
                mxDf=[];
                locmx=[];
                for i1=1:length(Loc_grd2mx1)
                    %             if vus1(Loc_grd2mx1(i1)+1,2)>0 | vus1(Loc_grd2mx1(i1)-1+1,2)>0 | vus1(Loc_grd2mx1(i1)+1+1,2)>0
                    if vus1(Loc_grd2mx1(i1)-1,2)>0
                        i2=i2+1;
                        mxDf(i2)=Grd2mx1(i1);
                        locmx(i2)=Loc_grd2mx1(i1);
                    end
                end
                
                i2=0;
                mnDf=[];
                locmn=[];
                for i1=1:length(Loc_grd2mn1)
                    %             if vus1(Loc_grd2mn1(i1)+1,2)>0 | vus1(Loc_grd2mn1(i1)-1+1,2)>0 | vus1(Loc_grd2mn1(i1)+1+1,2)>0
                    if vus1(Loc_grd2mn1(i1),2)>0
                        i2=i2+1;
                        mnDf(i2)=Grd2mn1(i1);
                        locmn(i2)=Loc_grd2mn1(i1);
                    end
                end
                
                LoC=sort([locmx locmn])';               % location of jumps at the SN point
                MxMn=[mxDf mnDf];
                jmpN=length(LoC);                       % number of jumps
                
                %**************************************************************************
                % Finding SN points from the unstable branch
                %**************************************************************************
                
                nonNanUS=(~isnan(vus1(:,2)));
%                 bfrPt=find(diff(sign(nonNanUS)))+1;      % SN points from the unstable branch
                bfrPt=find(diff(nonNanUS))+1;      % for 2015 verson of Matlab   
                chgSS=vss1(bfrPt,2)-vss1(bfrPt-1,2);
                bfrSN=[bfrPt chgSS];
                %         jmp=bfrSN(abs(bfrSN(:,2))>=4,2);
                %         bfr=bfrSN(abs(bfrSN(:,2))>=4,1);
                jmpLocSN=find(ismember(bfrSN(:,1),LoC));      % location of jump based on the ordering of SN points
                jmpSN=bfrSN(jmpLocSN,2);
                %**************************************************************************
                
                %         [max(peaks_ss) sum_diffss2 jmpN]
                
                %         if jmp > 0
                %             if length(loc_grd2mx1)~=length(loc_grd2mn1)
                %                 disp('yes')
                %                 continue
                %             end
                %         end
                
                % index for the parameter combination
                ic1=ic1+1;
                bistability_par(ic1,1)=kk;
                para_BS_pck(ic1,:)=parM(kk,:);
                
                % Plotting the bifurcations segregated according to their type
                
                %         fig2=figure(2);
                if(peaks_ss(1)>1 && peaks_ss(end)>1)          % irreversible switch on Left & right
                    if (sum_diffss2==0)
                        disp('LRbs')
                        LRbs=LRbs+1;
                        LRbs_par(LRbs,1)=kk;
                        para_LRbs(LRbs,:)=parM(kk,:);
                        phaseIndx(kp1,kp2)=116;
                    elseif (sum_diffss2==2)
                        if jmpN==0
                            disp('LRisis')
                            LRisis=LRisis+1;
                            LRisis_par(LRisis,1)=kk;
                            para_LRisis(LRisis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=117;
                        elseif jmpN==1
                            if jmpLocSN(1)==1
                                disp('LRinisis')
                                LRinisis=LRinisis+1;
                                LRinisis_par(LRinisis,1)=kk;
                                para_LRinisis(LRinisis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=118;
                            elseif jmpLocSN(1)==2
                                disp('LRisinis')
                                LRisinis=LRisinis+1;
                                LRisinis_par(LRisinis,1)=kk;
                                para_LRisinis(LRisinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=119;
                            end
                        elseif jmpN==2
                            disp('LRinisinis')
                            LRinisinis=LRinisinis+1;
                            LRinisinis_par(LRinisinis,1)=kk;
                            para_LRinisinis(LRinisinis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=120;
                        end
                    elseif (sum_diffss2==4)
                        if jmpN==0
                            disp('LRisisis')
                            LRisisis=LRisisis+1;
                            LRisisis_par(LRisisis,1)=kk;
                            para_LRisisis(LRisisis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=121;
                        elseif jmpN==1
                            if jmpLocSN(1)==1
                                disp('LRinisisis')
                                LRinisisis=LRinisisis+1;
                                LRinisisis_par(LRinisisis,1)=kk;
                                para_LRinisisis(LRinisisis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=122;
                            elseif jmpLocSN(1)==2
                                disp('LRisbsBis')
                                LRisbsBis=LRisbsBis+1;
                                LRisbsBis_par(LRisbsBis,1)=kk;
                                para_LRisbsBis(LRisbsBis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=123;
                            elseif jmpLocSN(1)==3
                                disp('LRisbsFis')
                                LRisbsFis=LRisbsFis+1;
                                LRisbsFis_par(LRisbsFis,1)=kk;
                                para_LRisbsFis(LRisbsFis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=124;
                            elseif jmpLocSN(1)==4
                                disp('LRisisinis')
                                LRisisinis=LRisisinis+1;
                                LRisisinis_par(LRisisinis,1)=kk;
                                para_LRisisinis(LRisisinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=125;
                            end
                        elseif jmpN==2
                            if jmpLocSN(1)==1 && jmpLocSN(2)==2
                                disp('LRinisbsBis')
                                LRinisbsBis=LRinisbsBis+1;
                                LRinisbsBis_par(LRinisbsBis,1)=kk;
                                para_LRinisbsBis(LRinisbsBis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=126;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3
                                disp('LRinisbsFis')
                                LRinisbsFis=LRinisbsFis+1;
                                LRinisbsFis_par(LRinisbsFis,1)=kk;
                                para_LRinisbsFis(LRinisbsFis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=127;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==4
                                disp('LRinisisinis')
                                LRinisisinis=LRinisisinis+1;
                                LRinisisinis_par(LRinisisinis,1)=kk;
                                para_LRinisisinis(LRinisisinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=128;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3
                                disp('LRisinisis')
                                LRisinisis=LRisinisis+1;
                                LRisinisis_par(LRisinisis,1)=kk;
                                para_LRisinisis(LRisinisis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=129;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==4
                                disp('LRisbsBinis')
                                LRisbsBinis=LRisbsBinis+1;
                                LRisbsBinis_par(LRisbsBinis,1)=kk;
                                para_LRisbsBinis(LRisbsBinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=130;
                            elseif jmpLocSN(1)==3 && jmpLocSN(2)==4
                                disp('LRisbsFinis')
                                LRisbsFinis=LRisbsFinis+1;
                                LRisbsFinis_par(LRisbsFinis,1)=kk;
                                para_LRisbsFinis(LRisbsFinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=131;
                            end
                        elseif jmpN==3
                            if jmpLocSN(1)==1 && jmpLocSN(2)==2 && jmpLocSN(3)==3
                                disp('LRinisinisis')
                                LRinisinisis=LRinisinisis+1;
                                LRinisinisis_par(LRinisinisis,1)=kk;
                                para_LRinisinisis(LRinisinisis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=132;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==2 && jmpLocSN(3)==4
                                disp('LRinisbsBinis')
                                LRinisbsBinis=LRinisbsBinis+1;
                                LRinisbsBinis_par(LRinisbsBinis,1)=kk;
                                para_LRinisbsBinis(LRinisbsBinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=133;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3 && jmpLocSN(3)==4
                                disp('LRinisbsFinis')
                                LRinisbsFinis=LRinisbsFinis+1;
                                LRinisbsFinis_par(LRinisbsFinis,1)=kk;
                                para_LRinisbsFinis(LRinisbsFinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=134;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3 && jmpLocSN(3)==4
                                disp('LRisinisinis')
                                LRisinisinis=LRisinisinis+1;
                                LRisinisinis_par(LRisinisinis,1)=kk;
                                para_LRisinisinis(LRisinisinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=135;
                            end
                        elseif jmpN==4
                            disp('LRinisinisinis')
                                LRinisinisinis=LRinisinisinis+1;
                                LRinisinisinis_par(LRinisinisinis,1)=kk;
                                para_LRinisinisinis(LRinisinisinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=136;
                        end
                    end
                end
                    
                if(peaks_ss(1)>1)          % irreversible switch on Left
                    if (sum_diffss2==1)
                        if jmpN==0
                            disp('LbsB')
                            LbsB=LbsB+1;
                            LbsB_par(LbsB,1)=kk;
                            para_LbsB(LbsB,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=60;
                        elseif jmpN==1
                            disp('LbsF')
                            LbsF=LbsF+1;
                            LbsF_par(LbsF,1)=kk;
                            para_LbsF(LbsF,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=61;
                        end
                    elseif (sum_diffss2==3)
                        if jmpN==0
                            disp('Ldlis')
                            Ldlis=Ldlis+1;
                            Ldlis_par(Ldlis,1)=kk;
                            para_Ldlis(Ldlis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=62;
                        elseif jmpN==1
                            if jmpLocSN(1)==1
                                disp('LbsFis')
                                LbsFis=LbsFis+1;
                                LbsFis_par(LbsFis,1)=kk;
                                para_LbsFis(LbsFis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=63;
                            elseif jmpLocSN(1)==2
                                disp('LbsBis')
                                LbsBis=LbsBis+1;
                                LbsBis_par(LbsBis,1)=kk;
                                para_LbsBis(LbsBis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=64;
                            elseif jmpLocSN(1)==3
                                disp('Linmus')
                                Linmus=Linmus+1;
                                Linmus_par(Linmus,1)=kk;
                                para_Linmus(Linmus,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=65;
                            end
                        elseif jmpN==2
                            if jmpLocSN(1)==1 && jmpLocSN(2)==3
                                disp('LdlbsF')
                                LdlbsF=LdlbsF+1;
                                LdlbsF_par(LdlbsF,1)=kk;
                                para_LdlbsF(LdlbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=66;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==2
                                disp('Lmus')
                                Lmus=Lmus+1;
                                Lmus_par(Lmus,1)=kk;
                                para_Lmus(Lmus,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=67;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3
                                disp('LbsBinis')
                                LbsBinis=LbsBinis+1;
                                LbsBinis_par(LbsBinis,1)=kk;
                                para_LbsBinis(LbsBinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=68;
                            end
                        elseif jmpN==3
                            disp('LbsFinis')
                            LbsFinis=LbsFinis+1;
                            LbsFinis_par(LbsFinis,1)=kk;
                            para_LbsFinis(LbsFinis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=69;
                        end
                    elseif (sum_diffss2==5)
                        if jmpN==0
                            disp('LdlisbsB')
                            LdlisbsB=LdlisbsB+1;
                            LdlisbsB_par(LdlisbsB,1)=kk;
                            para_LdlisbsB(LdlisbsB,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=70;
                        elseif jmpN==1
                            if jmpLocSN(1)==1
                                disp('LdlisbsF')
                                LdlisbsF=LdlisbsF+1;
                                LdlisbsF_par(LdlisbsF,1)=kk;
                                para_LdlisbsF(LdlisbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=71;
                            elseif jmpLocSN(1)==2
                                disp('LdlbsBisa')
                                LdlbsBisa=LdlbsBisa+1;
                                LdlbsBisa_par(LdlbsBisa,1)=kk;
                                para_LdlbsBisa(LdlbsBisa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=72;
                            elseif jmpLocSN(1)==3
                                disp('Linmusisa')
                                Linmusisa=Linmusisa+1;
                                Linmusisa_par(Linmusisa,1)=kk;
                                para_Linmusisa(Linmusisa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=73;
                            elseif jmpLocSN(1)==4
                                disp('LdlbsBism')
                                LdlbsBism=LdlbsBism+1;
                                LdlbsBism_par(LdlbsBism,1)=kk;
                                para_LdlbsBism(LdlbsBism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=74;
                            elseif jmpLocSN(1)==5
                                disp('Linmusism')
                                Linmusism=Linmusism+1;
                                Linmusism_par(Linmusism,1)=kk;
                                para_Linmusism(Linmusism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=75;
                            end
                        elseif jmpN==2
                            if jmpLocSN(1)==1 && jmpLocSN(2)==2
                                disp('Lmusisa')
                                Lmusisa=Lmusisa+1;
                                Lmusisa_par(Lmusisa,1)=kk;
                                para_Lmusisa(Lmusisa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=76;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3
                                disp('LdlbsFisa')
                                LdlbsFisa=LdlbsFisa+1;
                                LdlbsFisa_par(LdlbsFisa,1)=kk;
                                para_LdlbsFisa(LdlbsFisa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=77;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==4
                                disp('Lmusism')
                                Lmusism=Lmusism+1;
                                Lmusism_par(Lmusism,1)=kk;
                                para_Lmusism(Lmusism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=78;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==5
                                disp('LdlbsFism')
                                LdlbsFism=LdlbsFism+1;
                                LdlbsFism_par(LdlbsFism,1)=kk;
                                para_LdlbsFism(LdlbsFism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=79;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3
                                disp('LbsBisinism')
                                LbsBisinism=LbsBisinism+1;
                                LbsBisinism_par(LbsBisinism,1)=kk;
                                para_LbsBisinism(LbsBisinism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=80;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==4
                                disp('LtpbsB')
                                LtpbsB=LtpbsB+1;
                                LtpbsB_par(LtpbsB,1)=kk;
                                para_LtpbsB(LtpbsB,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=81;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==5
                                disp('LbsBinmus')
                                LbsBinmus=LbsBinmus+1;
                                LbsBinmus_par(LbsBinmus,1)=kk;
                                para_LbsBinmus(LbsBinmus,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=82;
                            elseif jmpLocSN(1)==3 && jmpLocSN(2)==4
                                disp('LmusbsB')
                                LmusbsB=LmusbsB+1;
                                LmusbsB_par(LmusbsB,1)=kk;
                                para_LmusbsB(LmusbsB,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=83;
                            elseif jmpLocSN(1)==3 && jmpLocSN(2)==5
                                disp('LinmusbsFa')
                                LinmusbsFa=LinmusbsFa+1;
                                LinmusbsFa_par(LinmusbsFa,1)=kk;
                                para_LinmusbsFa(LinmusbsFa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=83;
                            elseif jmpLocSN(1)==4 && jmpLocSN(2)==5
                                disp('LbsBisinisa')
                                LbsBisinisa=LbsBisinisa+1;
                                LbsBisinisa_par(LbsBisinisa,1)=kk;
                                para_LbsBisinisa(LbsBisinisa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=85;
                            end
                        end
                    end
                elseif (peaks_ss(end)>1)   % irreversible switch on Right
                    if (sum_diffss2==1)
                        if jmpN==0
                            disp('RbsF')
                            RbsF=RbsF+1;
                            RbsF_par(RbsF,1)=kk;
                            para_RbsF(RbsF,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=90;
                        elseif jmpN==1
                            disp('RbsB')
                            RbsB=RbsB+1;
                            RbsB_par(RbsB,1)=kk;
                            para_RbsB(RbsB,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=91;
                        end
                    elseif (sum_diffss2==3)
                        if jmpN==0
                            disp('Rdlis')
                            Rdlis=Rdlis+1;
                            Rdlis_par(Rdlis,1)=kk;
                            para_Rdlis(Rdlis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=92;
                        elseif jmpN==1
                            if jmpLocSN(1)==3
                                disp('RbsBis')
                                RbsBis=RbsBis+1;
                                RbsBis_par(RbsBis,1)=kk;
                                para_RbsBis(RbsBis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=93;
                            elseif jmpLocSN(1)==2
                                disp('RbsFis')
                                RbsFis=RbsFis+1;
                                RbsFis_par(RbsFis,1)=kk;
                                para_RbsFis(RbsFis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=94;
                            elseif jmpLocSN(1)==1
                                disp('Rinmus')
                                Rinmus=Rinmus+1;
                                Rinmus_par(Rinmus,1)=kk;
                                para_Rinmus(Rinmus,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=95;
                            end
                        elseif jmpN==2
                            if jmpLocSN(1)==1 && jmpLocSN(2)==2
                                disp('RbsFinis')
                                RbsFinis=RbsFinis+1;
                                RbsFinis_par(RbsFinis,1)=kk;
                                para_RbsFinis(RbsFinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=96;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3
                                disp('Rmus')
                                Rmus=Rmus+1;
                                Rmus_par(Rmus,1)=kk;
                                para_Rmus(Rmus,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=97;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3
                                disp('RbsBinis')
                                RbsBinis=RbsBinis+1;
                                RbsBinis_par(RbsBinis,1)=kk;
                                para_RbsBinis(RbsBinis,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=98;
                            end
                        elseif jmpN==3
                            disp('Rdlinis')
                            Rdlinis=Rdlinis+1;
                            Rdlinis_par(Rdlinis,1)=kk;
                            para_Rdlinis(Rdlinis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=99;
                        end
                    elseif (sum_diffss2==5)
                        if jmpN==0
                            disp('RdlisbsF')
                            RdlisbsF=RdlisbsF+1;
                            RdlisbsF_par(RdlisbsF,1)=kk;
                            para_RdlisbsF(RdlisbsF,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=100;
                        elseif jmpN==1
                            if jmpLocSN(1)==1
                                disp('Rinmusism')
                                Rinmusism=Rinmusism+1;
                                Rinmusism_par(Rinmusism,1)=kk;
                                para_Rinmusism(Rinmusism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=101;
                            elseif jmpLocSN(1)==2
                                disp('RdlbsFism')
                                RdlbsFism=RdlbsFism+1;
                                RdlbsFism_par(RdlbsFism,1)=kk;
                                para_RdlbsFism(RdlbsFism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=102;
                            elseif jmpLocSN(1)==3
                                disp('Rinmusisb')
                                Rinmusisb=Rinmusisb+1;
                                Rinmusisb_par(Rinmusisb,1)=kk;
                                para_Rinmusisb(Rinmusisb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=103;
                            elseif jmpLocSN(1)==4
                                disp('RdlbsFisb')
                                RdlbsFisb=RdlbsFisb+1;
                                RdlbsFisb_par(RdlbsFisb,1)=kk;
                                para_RdlbsFisb(RdlbsFisb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=104;
                            elseif jmpLocSN(1)==5
                                disp('RdlisbsB')
                                RdlisbsB=RdlisbsB+1;
                                RdlisbsB_par(RdlisbsB,1)=kk;
                                para_RdlisbsB(RdlisbsB,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=105;
                            end
                        elseif jmpN==2
                            if jmpLocSN(1)==1 && jmpLocSN(2)==2
                                disp('RinisisbsF')
                                RinisisbsF=RinisisbsF+1;
                                RinisisbsF_par(RinisisbsF,1)=kk;
                                para_RinisisbsF(RinisisbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=106;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3
                                disp('RinmusbsB')
                                RinmusbsB=RinmusbsB+1;
                                RinmusbsB_par(RinmusbsB,1)=kk;
                                para_RinmusbsB(RinmusbsB,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=107;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==4
                                disp('RinmusbsF')
                                RinmusbsF=RinmusbsF+1;
                                RinmusbsF_par(RinmusbsF,1)=kk;
                                para_RinmusbsF(RinmusbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=108;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==5
                                disp('RdlbsBism')
                                RdlbsBism=RdlbsBism+1;
                                RdlbsBism_par(RdlbsBism,1)=kk;
                                para_RdlbsBism(RdlbsBism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=109;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3
                                disp('RmusbsF')
                                RmusbsF=RmusbsF+1;
                                RmusbsF_par(RmusbsF,1)=kk;
                                para_RmusbsF(RmusbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=110;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==4
                                disp('RtpbsF')
                                RtpbsF=RtpbsF+1;
                                RtpbsF_par(RtpbsF,1)=kk;
                                para_RtpbsF(RtpbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=111;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==5
                                disp('Rmusism')
                                Rmusism=Rmusism+1;
                                Rmusism_par(Rmusism,1)=kk;
                                para_Rmusism(Rmusism,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=112;
                            elseif jmpLocSN(1)==3 && jmpLocSN(2)==4
                                disp('RisinisbsF')
                                RisinisbsF=RisinisbsF+1;
                                RisinisbsF_par(RisinisbsF,1)=kk;
                                para_RisinisbsF(RisinisbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=113;
                            elseif jmpLocSN(1)==3 && jmpLocSN(2)==5
                                disp('RdlbsBisb')
                                RdlbsBisb=RdlbsBisb+1;
                                RdlbsBisb_par(RdlbsBisb,1)=kk;
                                para_RdlbsBisb(RdlbsBisb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=114;
                            elseif jmpLocSN(1)==4 && jmpLocSN(2)==5
                                disp('Rmusisb')
                                Rmusisb=Rmusisb+1;
                                Rmusisb_par(Rmusisb,1)=kk;
                                para_Rmusisb(Rmusisb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=115;
                            end
                        end
                    end
                else
                    
                    if (max(peaks_ss)==2 && (sum_diffss2==2))
                        if jmpN==0
                            disp('Normal Isola')
                            iso=iso+1;
                            iso_par(iso,1)=kk;
                            para_iso(iso,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=1;
                            %                 subplot(3,4,1);
                            %                 title('Normal Isola');
                            %                 plot(vss(:,1),vss(:,2),'.','color','k','MarkerSize',4);
                            %                 hold on;
                            %                 plot(vss(:,1),vss(:,3),'.','color','k','MarkerSize',4);
                            %                 plot(vuss(:,1),vuss(:,2),'.','color','r','MarkerSize',4);
                        elseif jmpN==1
                            if jmpSN>0
                                disp('Bistable Forward')
                                bsF=bsF+1;
                                bsF_par(bsF,1)=kk;
                                para_bsF(bsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=2;
                            elseif jmpSN<0
                                disp('Bistable Backward')
                                bsB=bsB+1;
                                bsB_par(bsB,1)=kk;
                                para_bsB(bsB,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=3;
                            end
                        elseif jmpN==2
                            disp('Inverted Isola')
                            inis=inis+1;
                            inis_par(inis,1)=kk;
                            para_inis(inis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=4;
                        end
                    elseif (max(peaks_ss)==2 && (sum_diffss2==4))
                        if jmpN==0
                            disp('Dual Normal Isola')
                            dlis=dlis+1;
                            dlis_par(dlis,1)=kk;
                            para_dlis(dlis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=5;
                        elseif jmpN==1
                            if jmpLocSN==1
                                disp('Isola BS Bkwrd Aftr')
                                isbsBa=isbsBa+1;
                                isbsBa_par(isbsBa,1)=kk;
                                para_isbsBa(isbsBa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=6;
                            elseif jmpLocSN==2
                                disp('Isola BS Frwrd Aftr')
                                isbsFa=isbsFa+1;
                                isbsFa_par(isbsFa,1)=kk;
                                para_isbsFa(isbsFa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=7;
                            elseif jmpLocSN==3
                                disp('Isola BS Bkwrd Bfr')
                                isbsBb=isbsBb+1;
                                isbsBb_par(isbsBb,1)=kk;
                                para_isbsBb(isbsBb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=8;
                            elseif jmpLocSN==4
                                disp('Isola BS Frwrd Bfr')
                                isbsFb=isbsFb+1;
                                isbsFb_par(isbsFb,1)=kk;
                                para_isbsFb(isbsFb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=9;
                            end
                        elseif jmpN==2
                            
                            if jmpLocSN(1)==2 && jmpLocSN(2)==4
                                disp('Dual Bistable Frwrd')
                                dlbsF=dlbsF+1;
                                dlbsF_par(dlbsF,1)=kk;
                                para_dlbsF(dlbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=10;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3
                                disp('Dual Bistable Bkwrd')
                                dlbsB=dlbsB+1;
                                dlbsB_par(dlbsB,1)=kk;
                                para_dlbsB(dlbsB,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=11;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3
                                disp('Mushroom')
                                mus=mus+1;
                                mus_par(mus,1)=kk;
                                para_mus(mus,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=12;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==4
                                disp('Inverted Mushroom')
                                inmus=inmus+1;
                                inmus_par(inmus,1)=kk;
                                para_inmus(inmus,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=13;
                            end
                        elseif jmpN==3
                            if jmpLocSN(1)==1 && jmpLocSN(2)==2 && jmpLocSN(3)==4
                                disp('Inverted Isola BS Frwrd Bfr')
                                inisbsFb=inisbsFb+1;
                                inisbsFb_par(inisbsFb,1)=kk;
                                para_inisbsFb(inisbsFb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=14;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3 && jmpLocSN(3)==4
                                disp('Inverted Isola BS Frwrd Aftr')
                                inisbsFa=inisbsFa+1;
                                inisbsFa_par(inisbsFa,1)=kk;
                                para_inisbsFa(inisbsFa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=15;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==2 && jmpLocSN(3)==3
                                disp('Inverted Isola BS Bkwrd Bfr')
                                inisbsBb=inisbsBb+1;
                                inisbsBb_par(inisbsBb,1)=kk;
                                para_inisbsBb(inisbsBb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=16;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3 && jmpLocSN(3)==4
                                disp('Inverted Isola BS Bkwrd Aftr')
                                inisbsBa=inisbsBa+1;
                                inisbsBa_par(inisbsBa,1)=kk;
                                para_inisbsBa(inisbsBa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=17;
                            end
                        elseif jmpN==4
                            disp('Dual Inverted Isola')
                            dlinis=dlinis+1;
                            dlinis_par(dlinis,1)=kk;
                            para_dlinis(dlinis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=18;
                        end
                    elseif (max(peaks_ss)==2 && (sum_diffss2==6))
                        if jmpN==0
                            disp('Triple Isola')
                            tlis=tlis+1;
                            tlis_par(tlis,1)=kk;
                            para_tlis(tlis,:)=parM(kk,:);
                            phaseIndx(kp1,kp2)=19;
                        elseif jmpN==1
                            if jmpLocSN(1)==6
                                disp('Dual Isola BS Frwrd bb')
                                dlisbsFbb=dlisbsFbb+1;
                                dlisbsFbb_par(dlisbsFbb,1)=kk;
                                para_dlisbsFbb(dlisbsFbb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=20;
                            elseif jmpLocSN(1)==4
                                disp('Dual Isola BS Frwrd ba')
                                dlisbsFba=dlisbsFba+1;
                                dlisbsFba_par(dlisbsFba,1)=kk;
                                para_dlisbsFba(dlisbsFba,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=21;
                            elseif jmpLocSN(1)==2
                                disp('Dual Isola BS Frwrd aa')
                                dlisbsFaa=dlisbsFaa+1;
                                dlisbsFaa_par(dlisbsFaa,1)=kk;
                                para_dlisbsFaa(dlisbsFaa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=22;
                            elseif jmpLocSN(1)==5
                                disp('Dual Isola BS Bkwrd bb')
                                dlisbsBbb=dlisbsBbb+1;
                                dlisbsBbb_par(dlisbsBbb,1)=kk;
                                para_dlisbsBbb(dlisbsBbb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=23;
                            elseif jmpLocSN(1)==3
                                disp('Dual Isola BS Bkwrd ba')
                                dlisbsBba=dlisbsBba+1;
                                dlisbsBba_par(dlisbsBba,1)=kk;
                                para_dlisbsBba(dlisbsBba,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=24;
                            elseif jmpLocSN(1)==1
                                disp('Dual Isola BS Bkwrd aa')
                                dlisbsBaa=dlisbsBaa+1;
                                dlisbsBaa_par(dlisbsBaa,1)=kk;
                                para_dlisbsBaa(dlisbsBaa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=25;
                            end
                        elseif jmpN==2
                            if jmpSN(1)>0 && jmpSN(2)<0
                                if jmpLocSN(1)==4 && jmpLocSN(2)==5
                                    disp('Isola Mushroom Bfr')
                                    ismusb=ismusb+1;
                                    ismusb_par(ismusb,1)=kk;
                                    para_ismusb(ismusb,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=26;
                                elseif jmpLocSN(1)==2 && jmpLocSN(2)==5
                                    disp('Isola Mushroom Mdl')
                                    ismusm=ismusm+1;
                                    ismusm_par(ismusm,1)=kk;
                                    para_ismusm(ismusm,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=27;
                                elseif jmpLocSN(1)==2 && jmpLocSN(2)==3
                                    disp('Isola Mushroom Aftr')
                                    ismusa=ismusa+1;
                                    ismusa_par(ismusa,1)=kk;
                                    para_ismusa(ismusa,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=28;
                                end
                            elseif jmpSN(1)>0 && jmpSN(2)>0
                                if jmpLocSN(1)==4 && jmpLocSN(2)==6
                                    disp('Isola DBS Frwrd Bfr')
                                    isdlbsFb=isdlbsFb+1;
                                    isdlbsFb_par(isdlbsFb,1)=kk;
                                    para_isdlbsFb(isdlbsFb,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=29;
                                elseif jmpLocSN(1)==2 && jmpLocSN(2)==6
                                    disp('Isola DBS Frwrd Mdl')
                                    isdlbsFm=isdlbsFm+1;
                                    isdlbsFm_par(isdlbsFm,1)=kk;
                                    para_isdlbsFm(isdlbsFm,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=30;
                                elseif jmpLocSN(1)==2 && jmpLocSN(2)==4
                                    disp('Isola DBS Frwrd Aftr')
                                    isdlbsFa=isdlbsFa+1;
                                    isdlbsFa_par(isdlbsFa,1)=kk;
                                    para_isdlbsFa(isdlbsFa,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=31;
                                end
                            elseif jmpSN(1)<0 && jmpSN(2)<0
                                if jmpLocSN(1)==3 && jmpLocSN(2)==5
                                    disp('Isola DBS Frwrd Bfr')
                                    isdlbsBb=isdlbsBb+1;
                                    isdlbsBb_par(isdlbsBb,1)=kk;
                                    para_isdlbsBb(isdlbsBb,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=32;
                                elseif jmpLocSN(1)==1 && jmpLocSN(2)==5
                                    disp('Isola DBS Frwrd Mdl')
                                    isdlbsBm=isdlbsBm+1;
                                    isdlbsBm_par(isdlbsBm,1)=kk;
                                    para_isdlbsBm(isdlbsBm,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=33;
                                elseif jmpLocSN(1)==1 && jmpLocSN(2)==3
                                    disp('Isola DBS Frwrd Aftr')
                                    isdlbsBa=isdlbsBa+1;
                                    isdlbsBa_par(isdlbsBa,1)=kk;
                                    para_isdlbsBa(isdlbsBa,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=34;
                                end
                            elseif jmpSN(1)<0 && jmpSN(2)>0
                                if jmpLocSN(1)==3 && jmpLocSN(2)==6
                                    disp('Isola inverted mushroom Bfr')
                                    isinmusb=isinmusb+1;
                                    isinmusb_par(isinmusb,1)=kk;
                                    para_isinmusb(isinmusb,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=35;
                                elseif jmpLocSN(1)==1 && jmpLocSN(2)==6
                                    disp('Isola inverted mushroom Mdl')
                                    isinmusm=isinmusm+1;
                                    isinmusm_par(isinmusm,1)=kk;
                                    para_isinmusm(isinmusm,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=36;
                                elseif jmpLocSN(1)==1 && jmpLocSN(2)==4
                                    disp('Isola inverted mushroom Aftr')
                                    isinmusa=isinmusa+1;
                                    isinmusa_par(isinmusa,1)=kk;
                                    para_isinmusa(isinmusa,:)=parM(kk,:);
                                    phaseIndx(kp1,kp2)=37;
                                end
                            end
                        elseif jmpN==3
                            if jmpLocSN(1)==2 && jmpLocSN(2)==4 && jmpLocSN(3)==6
                                disp('Triple BS Frwrd')
                                tpbsF=tpbsF+1;
                                tpbsF_par(tpbsF,1)=kk;
                                para_tpbsF(tpbsF,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=38;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3 && jmpLocSN(3)==5
                                disp('Triple BS Bkwrd')
                                tpbsB=tpbsB+1;
                                tpbsB_par(tpbsB,1)=kk;
                                para_tpbsB(tpbsB,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=39;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3 && jmpLocSN(3)==6
                                disp('Mushroom BS Frwrd Aftr')
                                musbsFa=musbsFa+1;
                                musbsFa_par(musbsFa,1)=kk;
                                para_musbsFa(musbsFa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=40;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3 && jmpLocSN(3)==5
                                disp('Mushroom BS Bkwrd Aftr')
                                musbsBa=musbsBa+1;
                                musbsBa_par(musbsBa,1)=kk;
                                para_musbsBa(musbsBa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=41;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==4 && jmpLocSN(3)==5
                                disp('Mushroom BS Frwrd Bfr')
                                musbsFb=musbsFb+1;
                                musbsFb_par(musbsFb,1)=kk;
                                para_musbsFb(musbsFb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=42;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==4 && jmpLocSN(3)==5
                                disp('Inverted Mushroom BS Bkwrd Aftr')
                                inmusbsBa=inmusbsBa+1;
                                inmusbsBa_par(inmusbsBa,1)=kk;
                                para_inmusbsBa(inmusbsBa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=43;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==4 && jmpLocSN(3)==6
                                disp('Inverted Mushroom BS Frwrd Aftr')
                                inmusbsFa=inmusbsFa+1;
                                inmusbsFa_par(inmusbsFa,1)=kk;
                                para_inmusbsFa(inmusbsFa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=44;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3 && jmpLocSN(3)==6
                                disp('Inverted Mushroom BS Bkwrd Bfr')
                                inmusbsBb=inmusbsBb+1;
                                inmusbsBb_par(inmusbsBb,1)=kk;
                                para_inmusbsBb(inmusbsBb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=45;
                            end
                        elseif jmpN==4
                            if jmpLocSN(1)==1 && jmpLocSN(2)==2 && jmpLocSN(3)==3 && jmpLocSN(4)==6
                                disp('Inverted Isola Inverted Mushroom Bfr')
                                inisinmusb=inisinmusb+1;
                                inisinmusb_par(inisinmusb,1)=kk;
                                para_inisinmusb(inisinmusb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=46;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3 && jmpLocSN(3)==4 && jmpLocSN(4)==6
                                disp('Inverted Isola Inverted Mushroom Mdl')
                                inisinmusm=inisinmusm+1;
                                inisinmusm_par(inisinmusm,1)=kk;
                                para_inisinmusm(inisinmusm,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=47;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==4 && jmpLocSN(3)==5 && jmpLocSN(4)==6
                                disp('Inverted Isola Inverted Mushroom Aftr')
                                inisinmusa=inisinmusa+1;
                                inisinmusa_par(inisinmusa,1)=kk;
                                para_inisinmusa(inisinmusa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=48;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==2 && jmpLocSN(3)==4 && jmpLocSN(4)==5
                                disp('Inverted Isola Mushroom Bfr')
                                inismusb=inismusb+1;
                                inismusb_par(inismusb,1)=kk;
                                para_inismusb(inismusb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=49;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3 && jmpLocSN(3)==4 && jmpLocSN(4)==5
                                disp('Inverted Isola Mushroom Mdl')
                                inismusm=inismusm+1;
                                inismusm_par(inismusm,1)=kk;
                                para_inismusm(inismusm,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=50;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3 && jmpLocSN(3)==5 && jmpLocSN(4)==6
                                disp('Inverted Isola Mushroom Aftr')
                                inismusa=inismusa+1;
                                inismusa_par(inismusa,1)=kk;
                                para_inismusa(inismusa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=51;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==2 && jmpLocSN(3)==4 && jmpLocSN(4)==6
                                disp('Inverted Isola dual BS Frwrd Bfr')
                                inisdlbsFb=inisdlbsFb+1;
                                inisdlbsFb_par(inisdlbsFb,1)=kk;
                                para_inisdlbsFb(inisdlbsFb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=52;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==3 && jmpLocSN(3)==4 && jmpLocSN(4)==6
                                disp('Inverted Isola dual BS Frwrd Mdl')
                                inisdlbsFm=inisdlbsFm+1;
                                inisdlbsFm_par(inisdlbsFm,1)=kk;
                                para_inisdlbsFm(inisdlbsFm,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=53;
                            elseif jmpLocSN(1)==2 && jmpLocSN(2)==4 && jmpLocSN(3)==5 && jmpLocSN(4)==6
                                disp('Inverted Isola dual BS Frwrd Aftr')
                                inisdlbsFa=inisdlbsFa+1;
                                inisdlbsFa_par(inisdlbsFa,1)=kk;
                                para_inisdlbsFa(inisdlbsFa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=54;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==2 && jmpLocSN(3)==3 && jmpLocSN(4)==5
                                disp('Inverted Isola dual BS Bkwrd Bfr')
                                inisdlbsBb=inisdlbsBb+1;
                                inisdlbsBb_par(inisdlbsBb,1)=kk;
                                para_inisdlbsBb(inisdlbsBb,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=55;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3 && jmpLocSN(3)==4 && jmpLocSN(4)==5
                                disp('Inverted Isola dual BS Bkwrd Mdl')
                                inisdlbsBm=inisdlbsBm+1;
                                inisdlbsBm_par(inisdlbsBm,1)=kk;
                                para_inisdlbsBm(inisdlbsBm,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=56;
                            elseif jmpLocSN(1)==1 && jmpLocSN(2)==3 && jmpLocSN(3)==5 && jmpLocSN(4)==6
                                disp('Inverted Isola dual BS Bkwrd Aftr')
                                inisdlbsBa=inisdlbsBa+1;
                                inisdlbsBa_par(inisdlbsBa,1)=kk;
                                para_inisdlbsBa(inisdlbsBa,:)=parM(kk,:);
                                phaseIndx(kp1,kp2)=57;
                            end
                        end
                    end
                end
                
                %         set(fig2,'visible','off')
                
%                 clear stbss unsss vss vuss vss1 vssN vus1 Grd2mn Grd2mn1 Grd2mx ...
%                     Grd2mx1 Loc_grd2mn Loc_grd2mn1 Loc_grd2mx Loc_grd2mx1 ...
%                     mxDf mnDf locmn locmx LoC MxMn
                
            end
        end
    end
end

% set(fig1,'Visible','on');
% set(fig2,'Visible','on');

tot_iso=iso;
tot_bsF=bsF;
tot_bsB=bsB;
tot_inis=inis;
tot_dlis=dlis;
tot_isbsBa=isbsBa;
tot_isbsFa=isbsFa;
tot_isbsBb=isbsBb;
tot_isbsFb=isbsFb;
tot_dlbsF=dlbsF;
tot_dlbsB=dlbsB;
tot_mus=mus;
tot_inmus=inmus;
tot_inisbsFb=inisbsFb;
tot_inisbsFa=inisbsFa;
tot_inisbsBb=inisbsBb;
tot_inisbsBa=inisbsBa;
tot_dlinis=dlinis;
tot_tlis=tlis;
tot_dlisbsFbb=dlisbsFbb;
tot_dlisbsFba=dlisbsFba;
tot_dlisbsFaa=dlisbsFaa;
tot_dlisbsBbb=dlisbsBbb;
tot_dlisbsBba=dlisbsBba;
tot_dlisbsBaa=dlisbsBaa;
tot_ismusb=ismusb;
tot_ismusm=ismusm;
tot_ismusa=ismusa;
tot_isdlbsFb=isdlbsFb;
tot_isdlbsFm=isdlbsFm;
tot_isdlbsFa=isdlbsFa;
tot_isinmusb=isinmusb;
tot_isinmusm=isinmusm;
tot_isinmusa=isinmusa;
tot_tpbsF=tpbsF;
tot_tpbsB=tpbsB;
tot_musbsFa=musbsFa;
tot_musbsBa=musbsBa;
tot_musbsFb=musbsFb;
tot_inmusbsBa=inmusbsBa;
tot_inmusbsFa=inmusbsFa;
tot_inmusbsBb=inmusbsBb;
tot_inisinmusb=inisinmusb;
tot_inisinmusm=inisinmusm;
tot_inisinmusa=inisinmusa;
tot_inismusb=inismusb;
tot_inismusm=inismusm;
tot_inismusa=inismusa;
tot_inisdlbsFb=inisdlbsFb;
tot_inisdlbsFm=inisdlbsFm;
tot_inisdlbsFa=inisdlbsFa;
tot_inisdlbsBb=inisdlbsBb;
tot_inisdlbsBm=inisdlbsBm;
tot_inisdlbsBa=inisdlbsBa;
tot_LbsB=LbsB;
tot_LbsF=LbsF;
tot_Ldlis=Ldlis;
tot_LbsFis=LbsFis;
tot_LbsBis=LbsBis;
tot_Linmus=Linmus;
tot_LdlbsF=LdlbsF;
tot_Lmus=Lmus;
tot_LbsBinis=LbsBinis;
tot_LbsFinis=LbsFinis;
tot_RbsF=RbsF;
tot_RbsB=RbsB;
tot_Rdlis=Rdlis;
tot_RbsBis=RbsBis;
tot_RbsFis=RbsFis;
tot_Rinmus=Rinmus;
tot_RbsFinis=RbsFinis;
tot_Rmus=Rmus;
tot_RbsBinis=RbsBinis;
tot_Rdlinis=Rdlinis;
tot_LdlisbsB=LdlisbsB;
tot_LdlisbsF=LdlisbsF;
tot_LdlbsBisa=LdlbsBisa;
tot_Linmusisa=Linmusisa;
tot_LdlbsBism=LdlbsBism;
tot_Linmusism=Linmusism;
tot_Lmusisa=Lmusisa;
tot_LdlbsFisa=LdlbsFisa;
tot_Lmusism=Lmusism;
tot_LdlbsFism=LdlbsFism;
tot_LbsBisinism=LbsBisinism;
tot_LtpbsB=LtpbsB;
tot_LbsBinmus=LbsBinmus;
tot_LmusbsB=LmusbsB;
tot_LinmusbsFa=LinmusbsFa;
tot_LbsBisinisa=LbsBisinisa;
tot_RdlisbsF=RdlisbsF;
tot_Rinmusism=Rinmusism;
tot_RdlbsFism=RdlbsFism;
tot_Rinmusisb=Rinmusisb;
tot_RdlbsFisb=RdlbsFisb;
tot_RdlisbsB=RdlisbsB;
tot_RinisisbsF=RinisisbsF;
tot_RinmusbsB=RinmusbsB;
tot_RinmusbsF=RinmusbsF;
tot_RdlbsBism=RdlbsBism;
tot_RmusbsF=RmusbsF;
tot_RtpbsF=RtpbsF;
tot_Rmusism=Rmusism;
tot_RisinisbsF=RisinisbsF;
tot_RdlbsBisb=RdlbsBisb;
tot_Rmusisb=Rmusisb;
tot_LRbs=LRbs;
tot_LRisis=LRisis;
tot_LRinisis=LRinisis;
tot_LRisinis=LRisinis;
tot_LRinisinis=LRinisinis;
tot_LRisisis=LRisisis;
tot_LRinisisis=LRinisisis;
tot_LRisbsBis=LRisbsBis;
tot_LRisbsFis=LRisbsFis;
tot_LRisisinis=LRisisinis;
tot_LRinisbsBis=LRinisbsBis;
tot_LRinisbsFis=LRinisbsFis;
tot_LRinisisinis=LRinisisinis;
tot_LRisinisis=LRisinisis;
tot_LRisbsBinis=LRisbsBinis;
tot_LRisbsFinis=LRisbsFinis;
tot_LRinisinisis=LRinisinisis;
tot_LRinisbsBinis=LRinisbsBinis;
tot_LRinisbsFinis=LRinisbsFinis;
tot_LRisinisinis=LRisinisinis;
tot_LRinisinisinis=LRinisinisinis;

% ppMISA_OR_dataR.total_bistability=[tot_bis tot_iso tot_bsF tot_bsB tot_inis ...
%     tot_dlis tot_isbsBa tot_isbsFa tot_isbsBb tot_isbsFb tot_inis tot_dlis ...
%     tot_dlbsF tot_dlbsB tot_mus=mus tot_inmus=inmus tot_inisbsFb ...
%     tot_inisbsFa tot_inisbsBb tot_inisbsBa tot_dlinis tot_tlis ...
%     tot_dlisbsFbb tot_dlisbsFba tot_dlisbsFaa tot_dlisbsBbb ...
%     tot_dlisbsBba tot_dlisbsBaa tot_ismusb tot_ismusm tot_ismusa ...
%     tot_isdlbsFb tot_isdlbsFm tot_isdlbsFa tot_isdlbsBb tot_isdlbsBm tot_isdlbsBa ...
%     tot_isinmusb tot_isinmusm tot_isinmusa tot_tpbsF tot_tpbsB ...
%     tot_musbsFa tot_musbsBa tot_musbsFb tot_inmusbsBa tot_inmusbsFa ...
%     tot_inmusbsBb tot_inisinmusb tot_inisinmusm tot_inisinmusb ...
%     tot_inismusb tot_inismusm tot_inismusa tot_inisdlbsFb ...
%     tot_inisdlbsFm tot_inisdlbsFa tot_inisdlbsBb tot_inisdlbsBm ...
%     tot_inisdlbsFa tot_LbsB tot_Ldlis tot_LbsFis tot_LbsBis ...
%     tot_Linmus tot_LdlbsF Lmus LbsBinis LbsFinis tot_RbsF tot_RbsB ...
%     tot_Rdlis tot_RbsBis tot_RbsFis tot_Rinmus tot_RbsFinis ...
%     tot_Rmus tot_RbsBinis tot_Rdlinis tot_LdlisbsB tot_LdlisbsF ...
%     tot_LdlbsBisa tot_Linmusisa tot_LdlbsBism tot_Lmusisa tot_LdlbsFisa ...
%     tot_Lmusism tot_LdlbsFism tot_LbsBisinism tot_LtpbsB tot_LbsBinmus ...
%     tot_LmusbsB tot_LinmusbsFa tot_LbsBisinisa tot_RdlisbsF tot_Rinmusism ...
%     tot_RdlbsFism tot_Rinmusisb tot_RdlbsFisb tot_RdlisbsB tot_RinisisbsF ...
%     tot_RinmusbsB tot_RinmusbsF tot_RdlbsBism tot_RmusbsF tot_RtpbsF ...
%     tot_Rmusism tot_RisinisbsF tot_RdlbsBisb tot_Rmusisb tot_LRbs ...
%     tot_LRisis tot_LRinisis tot_LRisinis tot_LRinisinis tot_LRisisis ...
%     tot_LRinisisis tot_LRisbsBis tot_LRisbsFis tot_LRisisinis ...
%     tot_LRinisbsBis tot_LRinisbsFis tot_LRinisisinis tot_LRisinisis ...
%     tot_LRisbsBinis tot_LRisbsFinis tot_LRinisinisis tot_LRinisbsBinis ...
%     tot_LRinisbsFinis tot_LRisinisinis tot_LRinisinisinis];
 
% total_seg=sum(ppMISA_OR_dataR.total_bistability)
% discrepency=total_seg-length(bistability_par)
%
% ppMISA_OR_dataR.data=parM;
% ppMISA_OR_dataR.para_BS_all=parM;
% ppMISA_OR_dataR.para_BS_pck=parM;

% ppMISA_OR_dataR.para_iso=para_iso;
% ppMISA_OR_dataR.para_bsF=para_bsF;
% ppMISA_OR_dataR.para_bsB=para_bsB;
% ppMISA_OR_dataR.para_inis=para_inis;
% ppMISA_OR_dataR.para_dlis=para_dlis;
% ppMISA_OR_dataR.para_isbsBa=para_isbsBa;
% ppMISA_OR_dataR.para_isbsFa=para_isbsFa;
% ppMISA_OR_dataR.para_isbsBb=para_isbsBb;
% ppMISA_OR_dataR.para_isbsFb=para_isbsFb;
% ppMISA_OR_dataR.para_dlbsF=para_dlbsF;
% ppMISA_OR_dataR.para_dlbsB=para_dlbsB;
% ppMISA_OR_dataR.para_mus=para_mus;
% ppMISA_OR_dataR.para_inmus=para_inmus;
% ppMISA_OR_dataR.para_inisbsFb=para_inisbsFb;
% ppMISA_OR_dataR.para_inisbsFa=para_inisbsFa;
% ppMISA_OR_dataR.para_inisbsBb=para_inisbsBb;
% ppMISA_OR_dataR.para_inisbsBa=para_inisbsBa;
% ppMISA_OR_dataR.para_dlinis=para_dlinis;
% ppMISA_OR_dataR.para_tlis=para_tlis;
% ppMISA_OR_dataR.para_dlisbsFbb=para_dlisbsFbb;
% ppMISA_OR_dataR.para_dlisbsFba=para_dlisbsFba;
% ppMISA_OR_dataR.para_dlisbsFaa=para_dlisbsFaa;
% ppMISA_OR_dataR.para_dlisbsBbb=para_dlisbsBbb;
% ppMISA_OR_dataR.para_dlisbsBba=para_dlisbsBba;
% ppMISA_OR_dataR.para_dlisbsBaa=para_dlisbsBaa;
% ppMISA_OR_dataR.para_ismusb=para_ismusb;
% ppMISA_OR_dataR.para_ismusm=para_ismusm;
% ppMISA_OR_dataR.para_ismusa=para_ismusa;
% ppMISA_OR_dataR.para_isdlbsFb=para_isdlbsFb;
% ppMISA_OR_dataR.para_isdlbsFm=para_isdlbsFm;
% ppMISA_OR_dataR.para_isdlbsFa=para_isdlbsFa;
% ppMISA_OR_dataR.para_isdlbsBb=para_isdlbsBb;
% ppMISA_OR_dataR.para_isdlbsBm=para_isdlbsBm;
% ppMISA_OR_dataR.para_isdlbsBa=para_isdlbsBa;
% ppMISA_OR_dataR.para_isinmusb=para_isinmusb;
% ppMISA_OR_dataR.para_isinmusm=para_isinmusm;
% ppMISA_OR_dataR.para_isinmusa=para_isinmusa;
% ppMISA_OR_dataR.para_tpbsF=para_tpbsF;
% ppMISA_OR_dataR.para_tpbsB=para_tpbsB;
% ppMISA_OR_dataR.para_musbsFa=para_musbsFa;
% ppMISA_OR_dataR.para_musbsBa=para_musbsBa;
% ppMISA_OR_dataR.para_musbsFb=para_musbsFb;
% ppMISA_OR_dataR.para_inmusbsBa=para_inmusbsBa;
% ppMISA_OR_dataR.para_inmusbsFa=para_inmusbsFa;
% ppMISA_OR_dataR.para_inmusbsBb=para_inmusbsBb;
% ppMISA_OR_dataR.para_inisinmusb=para_inisinmusb;
% ppMISA_OR_dataR.para_inisinmusm=para_inisinmusm;
% ppMISA_OR_dataR.para_inisinmusa=para_inisinmusa;
% ppMISA_OR_dataR.para_inismusb=para_inismusb;
% ppMISA_OR_dataR.para_inismusm=para_inismusm;
% ppMISA_OR_dataR.para_inismusa=para_inismusa;
% ppMISA_OR_dataR.para_inisdlbsFb=para_inisdlbsFb;
% ppMISA_OR_dataR.para_inisdlbsFm=para_inisdlbsFm;
% ppMISA_OR_dataR.para_inisdlbsFa=para_inisdlbsFa;
% ppMISA_OR_dataR.para_inisdlbsBb=para_inisdlbsBb;
% ppMISA_OR_dataR.para_inisdlbsBm=para_inisdlbsBm;
% ppMISA_OR_dataR.para_inisdlbsBa=para_inisdlbsBa;
% ppMISA_OR_dataR.para_LbsB=para_LbsB;
% ppMISA_OR_dataR.para_LbsF=para_LbsF;
% ppMISA_OR_dataR.para_Ldlis=para_Ldlis;
% ppMISA_OR_dataR.para_LbsFis=para_LbsFis;
% ppMISA_OR_dataR.para_LbsBis=para_LbsBis;
% ppMISA_OR_dataR.para_Linmus=para_Linmus;
% ppMISA_OR_dataR.para_LdlbsF=para_LdlbsF;
% ppMISA_OR_dataR.para_Lmus=para_Lmus;
% ppMISA_OR_dataR.para_LbsBinis=para_LbsBinis;
% ppMISA_OR_dataR.para_LbsFinis=para_LbsFinis;
% ppMISA_OR_dataR.para_RbsF=para_RbsF;
% ppMISA_OR_dataR.para_RbsB=para_RbsB;
% ppMISA_OR_dataR.para_Rdlis=para_Rdlis;
% ppMISA_OR_dataR.para_RbsBis=para_RbsBis;
% ppMISA_OR_dataR.para_LbsFis=para_LbsFis;
% ppMISA_OR_dataR.para_Rinmus=para_Rinmus;
% ppMISA_OR_dataR.para_RbsFinis=para_RbsFinis;
% ppMISA_OR_dataR.para_Rmus=para_Rmus;
% ppMISA_OR_dataR.para_RbsBinis=para_RbsBinis;
% ppMISA_OR_dataR.para_Rdlinis=para_Rdlinis;
% ppMISA_OR_dataR.para_LdlisbsB=para_LdlisbsB;
% ppMISA_OR_dataR.para_LdlisbsF=para_LdlisbsF;
% ppMISA_OR_dataR.para_LdlbsBisa=para_LdlbsBisa;
% ppMISA_OR_dataR.para_Linmusisa=para_Linmusisa;
% ppMISA_OR_dataR.para_LdlbsBism=para_LdlbsBism;
% ppMISA_OR_dataR.para_Linmusism=para_Linmusism;
% ppMISA_OR_dataR.para_Lmusisa=para_Lmusisa;
% ppMISA_OR_dataR.para_LdlbsFisa=para_LdlbsFisa;
% ppMISA_OR_dataR.para_Lmusism=para_Lmusism;
% ppMISA_OR_dataR.para_LdlbsFism=para_LdlbsFism;
% ppMISA_OR_dataR.para_LbsBisinism=para_LbsBisinism;
% ppMISA_OR_dataR.para_LtpbsB=para_LtpbsB;
% ppMISA_OR_dataR.para_LbsBinmus=para_LbsBinmus;
% ppMISA_OR_dataR.para_LmusbsB=para_LmusbsB;
% ppMISA_OR_dataR.para_LinmusbsFa=para_LinmusbsFa;
% ppMISA_OR_dataR.para_LbsBisinisa=para_LbsBisinisa;
% ppMISA_OR_dataR.para_RdlisbsF=para_RdlisbsF;
% ppMISA_OR_dataR.para_Rinmusism=para_Rinmusism;
% ppMISA_OR_dataR.para_RdlbsFism=para_RdlbsFism;
% ppMISA_OR_dataR.para_Rinmusisb=para_Rinmusisb;
% ppMISA_OR_dataR.para_RdlbsFisb=para_RdlbsFisb;
% ppMISA_OR_dataR.para_RdlisbsB=para_RdlisbsB;
% ppMISA_OR_dataR.para_RinisisbsF=para_RinisisbsF;
% ppMISA_OR_dataR.para_RinmusbsB=para_RinmusbsB;
% ppMISA_OR_dataR.para_RinmusbsF=para_RinmusbsF;
% ppMISA_OR_dataR.para_RdlbsBism=para_RdlbsBism;
% ppMISA_OR_dataR.para_RmusbsF=para_RmusbsF;
% ppMISA_OR_dataR.para_RtpbsF=para_RtpbsF;
% ppMISA_OR_dataR.para_Rmusism=para_Rmusism;
% ppMISA_OR_dataR.para_RisinisbsF=para_RisinisbsF;
% ppMISA_OR_dataR.para_RdlbsBisb=para_RdlbsBisb;
% ppMISA_OR_dataR.para_Rmusisb=para_Rmusisb;
% ppMISA_OR_dataR.para_LRbs=para_LRbs;
% ppMISA_OR_dataR.para_LRisis=para_LRisis;
% ppMISA_OR_dataR.para_LRinisis=para_LRinisis;
% ppMISA_OR_dataR.para_LRisinis=para_LRisinis;
% ppMISA_OR_dataR.para_LRinisinis=para_LRinisinis;
% ppMISA_OR_dataR.para_LRisisis=para_LRisisis;
% ppMISA_OR_dataR.para_LRinisisis=para_LRinisisis;
% ppMISA_OR_dataR.para_LRisbsBis=para_LRisbsBis;
% ppMISA_OR_dataR.para_LRisbsFis=para_LRisbsFis;
% ppMISA_OR_dataR.para_LRisisinis=para_LRisisinis;
% ppMISA_OR_dataR.para_LRinisbsBis=para_LRinisbsBis;
% ppMISA_OR_dataR.para_LRinisbsFis=para_LRinisbsFis;
% ppMISA_OR_dataR.para_LRinisisinis=para_LRinisisinis;
% ppMISA_OR_dataR.para_LRisinisis=para_LRisinisis;
% ppMISA_OR_dataR.para_LRisbsBinis=para_LRisbsBinis;
% ppMISA_OR_dataR.para_LRisbsFinis=para_LRisbsFinis;
% ppMISA_OR_dataR.para_LRinisinisis=para_LRinisinisis;
% ppMISA_OR_dataR.para_LRinisbsBinis=para_LRinisbsBinis;
% ppMISA_OR_dataR.para_LRinisbsFinis=para_LRinisbsFinis;
% ppMISA_OR_dataR.para_LRisbsFinis=para_LRisbsFinis;
% ppMISA_OR_dataR.para_LRinisinisinis=para_LRinisinisinis;

% ppMISA_OR_dataR.indx_iso=iso_par;
% ppMISA_OR_dataR.indx_bsF=bsF_par;
% ppMISA_OR_dataR.indx_bsB=bsB_par;
% ppMISA_OR_dataR.indx_inis=inis_par;
% ppMISA_OR_dataR.indx_dlis=dlis_par;
% ppMISA_OR_dataR.indx_isbsBa=isbsBa_par;
% ppMISA_OR_dataR.indx_isbsFa=isbsFa_par;
% ppMISA_OR_dataR.indx_isbsBb=isbsBb_par;
% ppMISA_OR_dataR.indx_isbsFb=isbsFb_par;
% ppMISA_OR_dataR.indx_dlbsF=dlbsF_par;
% ppMISA_OR_dataR.indx_dlbsB=dlbsB_par;
% ppMISA_OR_dataR.indx_mus=mus_par;
% ppMISA_OR_dataR.indx_inmus=inmus_par;
% ppMISA_OR_dataR.indx_inisbsFb=inisbsFb_par;
% ppMISA_OR_dataR.indx_inisbsFa=inisbsFa_par;
% ppMISA_OR_dataR.indx_inisbsBb=inisbsBb_par;
% ppMISA_OR_dataR.indx_inisbsBa=inisbsBa_par;
% ppMISA_OR_dataR.indx_dlinis=dlinis_par;
% ppMISA_OR_dataR.indx_tlis=tlis_par;
% ppMISA_OR_dataR.indx_dlisbsFbb=dlisbsFbb_par;
% ppMISA_OR_dataR.indx_dlisbsFba=dlisbsFba_par;
% ppMISA_OR_dataR.indx_dlisbsFaa=dlisbsFaa_par;
% ppMISA_OR_dataR.indx_dlisbsBbb=dlisbsBbb_par;
% ppMISA_OR_dataR.indx_dlisbsBba=dlisbsBba_par;
% ppMISA_OR_dataR.indx_dlisbsBaa=dlisbsBaa_par;
% ppMISA_OR_dataR.indx_ismusb=ismusb_par;
% ppMISA_OR_dataR.indx_ismusm=ismusm_par;
% ppMISA_OR_dataR.indx_ismusa=ismusa_par;
% ppMISA_OR_dataR.indx_isdlbsFb=isdlbsFb_par;
% ppMISA_OR_dataR.indx_isdlbsFm=isdlbsFm_par;
% ppMISA_OR_dataR.indx_isdlbsFa=isdlbsFa_par;
% ppMISA_OR_dataR.indx_isdlbsBb=isdlbsBb_par;
% ppMISA_OR_dataR.indx_isdlbsBm=isdlbsBm_par;
% ppMISA_OR_dataR.indx_isdlbsBa=isdlbsBa_par;
% ppMISA_OR_dataR.indx_isinmusb=isinmusb_par;
% ppMISA_OR_dataR.indx_isinmusm=isinmusm_par;
% ppMISA_OR_dataR.indx_isinmusa=isinmusa_par;
% ppMISA_OR_dataR.indx_tpbsF=tpbsF_par;
% ppMISA_OR_dataR.indx_tpbsB=tpbsB_par;
% ppMISA_OR_dataR.indx_musbsFa=musbsFa_par;
% ppMISA_OR_dataR.indx_musbsBa=musbsBa_par;
% ppMISA_OR_dataR.indx_musbsFb=musbsFb_par;
% ppMISA_OR_dataR.indx_inmusbsBa=inmusbsBa_par;
% ppMISA_OR_dataR.indx_inmusbsFa=inmusbsFa_par;
% ppMISA_OR_dataR.indx_inmusbsBb=inmusbsBb_par;
% ppMISA_OR_dataR.indx_inisinmusb=inisinmusb_par;
% ppMISA_OR_dataR.indx_inisinmusm=inisinmusm_par;
% ppMISA_OR_dataR.indx_inisinmusa=inisinmusa_par;
% ppMISA_OR_dataR.indx_inismusb=inismusb_par;
% ppMISA_OR_dataR.indx_inismusm=inismusm_par;
% ppMISA_OR_dataR.indx_inismusa=inismusa_par;
% ppMISA_OR_dataR.indx_inisdlbsFb=inisdlbsFb_par;
% ppMISA_OR_dataR.indx_inisdlbsFm=inisdlbsFm_par;
% ppMISA_OR_dataR.indx_inisdlbsFa=inisdlbsFa_par;
% ppMISA_OR_dataR.indx_inisdlbsBb=inisdlbsBb_par;
% ppMISA_OR_dataR.indx_inisdlbsBm=inisdlbsBm_par;
% ppMISA_OR_dataR.indx_inisdlbsBa=inisdlbsBa_par;
% ppMISA_OR_dataR.indx_LbsB=LbsB_par;
% ppMISA_OR_dataR.indx_LbsF=LbsF_par;
% ppMISA_OR_dataR.indx_Ldlis=Ldlis_par;
% ppMISA_OR_dataR.indx_LbsFis=LbsFis_par;
% ppMISA_OR_dataR.indx_LbsBis=LbsBis_par;
% ppMISA_OR_dataR.indx_Linmus=Linmus_par;
% ppMISA_OR_dataR.indx_LdlbsF=LdlbsF_par;
% ppMISA_OR_dataR.indx_Lmus=Lmus_par;
% ppMISA_OR_dataR.indx_LbsBinis=LbsBinis_par;
% ppMISA_OR_dataR.indx_LbsFinis=LbsFinis_par;
% ppMISA_OR_dataR.indx_RbsF=RbsF_par;
% ppMISA_OR_dataR.indx_RbsB=RbsB_par;
% ppMISA_OR_dataR.indx_Rdlis=Rdlis_par;
% ppMISA_OR_dataR.indx_RbsBis=RbsBis_par;
% ppMISA_OR_dataR.indx_LbsFis=LbsFis_par;
% ppMISA_OR_dataR.indx_Rinmus=Rinmus_par;
% ppMISA_OR_dataR.indx_RbsFinis=RbsFinis_par;
% ppMISA_OR_dataR.indx_Rmus=Rmus_par;
% ppMISA_OR_dataR.indx_RbsBinis=RbsBinis_par;
% ppMISA_OR_dataR.indx_Rdlinis=Rdlinis_par;
% ppMISA_OR_dataR.indx_LdlisbsB=LdlisbsB_par;
% ppMISA_OR_dataR.indx_LdlisbsF=LdlisbsF_par;
% ppMISA_OR_dataR.indx_LdlbsBisa=LdlbsBisa_par;
% ppMISA_OR_dataR.indx_Linmusisa=Linmusisa_par;
% ppMISA_OR_dataR.indx_LdlbsBism=LdlbsBism_par;
% ppMISA_OR_dataR.indx_Linmusism=Linmusism_par;
% ppMISA_OR_dataR.indx_Lmusisa=Lmusisa_par;
% ppMISA_OR_dataR.indx_LdlbsFisa=LdlbsFisa_par;
% ppMISA_OR_dataR.indx_Lmusism=Lmusism_par;
% ppMISA_OR_dataR.indx_LdlbsFism=LdlbsFism_par;
% ppMISA_OR_dataR.indx_LbsBisinism=LbsBisinism_par;
% ppMISA_OR_dataR.indx_LtpbsB=LtpbsB_par;
% ppMISA_OR_dataR.indx_LbsBinmus=LbsBinmus_par;
% ppMISA_OR_dataR.indx_LmusbsB=LmusbsB_par;
% ppMISA_OR_dataR.indx_LinmusbsFa=LinmusbsFa_par;
% ppMISA_OR_dataR.indx_LbsBisinisa=LbsBisinisa_par;
% ppMISA_OR_dataR.indx_RdlisbsF=RdlisbsF_par;
% ppMISA_OR_dataR.indx_Rinmusism=Rinmusism_par;
% ppMISA_OR_dataR.indx_RdlbsFism=RdlbsFism_par;
% ppMISA_OR_dataR.indx_Rinmusisb=Rinmusisb_par;
% ppMISA_OR_dataR.indx_RdlbsFisb=RdlbsFisb_par;
% ppMISA_OR_dataR.indx_RdlisbsB=RdlisbsB_par;
% ppMISA_OR_dataR.indx_RinisisbsF=RinisisbsF_par;
% ppMISA_OR_dataR.indx_RinmusbsB=RinmusbsB_par;
% ppMISA_OR_dataR.indx_RinmusbsF=RinmusbsF_par;
% ppMISA_OR_dataR.indx_RdlbsBism=RdlbsBism_par;
% ppMISA_OR_dataR.indx_RmusbsF=RmusbsF_par;
% ppMISA_OR_dataR.indx_RtpbsF=RtpbsF_par;
% ppMISA_OR_dataR.indx_Rmusism=Rmusism_par;
% ppMISA_OR_dataR.indx_RisinisbsF=RisinisbsF_par;
% ppMISA_OR_dataR.indx_RdlbsBisb=RdlbsBisb_par;
% ppMISA_OR_dataR.indx_Rmusisb=Rmusisb_par;
% ppMISA_OR_dataR.indx_LRbs=LRbs_par;
% ppMISA_OR_dataR.indx_LRisis=LRisis_par;
% ppMISA_OR_dataR.indx_LRinisis=LRinisis_par;
% ppMISA_OR_dataR.indx_LRisinis=LRisinis_par;
% ppMISA_OR_dataR.indx_LRinisinis=LRinisinis_par;
% ppMISA_OR_dataR.indx_LRisisis=LRisisis_par;
% ppMISA_OR_dataR.indx_LRinisisis=LRinisisis_par;
% ppMISA_OR_dataR.indx_LRisbsBis=LRisbsBis_par;
% ppMISA_OR_dataR.indx_LRisbsFis=LRisbsFis_par;
% ppMISA_OR_dataR.indx_LRisisinis=LRisisinis_par;
% ppMISA_OR_dataR.indx_LRinisbsBis=LRinisbsBis_par;
% ppMISA_OR_dataR.indx_LRinisbsFis=LRinisbsFis_par;
% ppMISA_OR_dataR.indx_LRinisisinis=LRinisisinis_par;
% ppMISA_OR_dataR.indx_LRisinisis=LRisinisis_par;
% ppMISA_OR_dataR.indx_LRisbsBinis=LRisbsBinis_par;
% ppMISA_OR_dataR.indx_LRisbsFinis=LRisbsFinis_par;
% ppMISA_OR_dataR.indx_LRinisinisis=LRinisinisis_par;
% ppMISA_OR_dataR.indx_LRinisbsBinis=LRinisbsBinis_par;
% ppMISA_OR_dataR.indx_LRinisbsFinis=LRinisbsFinis_par;
% ppMISA_OR_dataR.indx_LRisbsFinis=LRisbsFinis_par;
% ppMISA_OR_dataR.indx_LRinisinisinis=LRinisinisinis_par;


% save ppMISA_OR_dataRRN_R1.mat ppMISA_OR_dataR

save phaseDia_ppMISA_Iso_r1p12_OR_r1.mat phaseIndx Val_gab Val_gba
savefig(fig1,'phaseDia_ppMISA_Iso_r1p12_OR_r1')
saveas(fig1,'phaseDia_ppMISA_Iso_r1p12_OR_r1','png')

toc

%quit

