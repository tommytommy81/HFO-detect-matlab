%% detectorUMC stage 3

% =========================================================================
% *** Function UMC_IIIstage_detectorUMC
% ***
% *** automatic time-frequency algorithm for detection of HFOs
% *** for more details refer to the publication
% *** Fedele et al. ??? "Automatically detected residual fast ripples" 
% ***
% ***
% ***----------------------------------------------------------------------
% *** Analysis (1 channel at the time):
% *** 0. Preprocessing
%        Filter data in the range [hp lp]
%        Calculate Hilbert envelope of the band passed signal
% *** 1. Stage I
%        Calculate threshold according to entropy based method
%        Detection of Events of Interest
%        Merge EoIs with inter-event-interval less than 10 ms into one EoI
% *** 2. Stage II
%        Compute TF around the EOI and check for contribution in HF
%        Collect statistic of each event

 
% -----INPUT---------------------------------------------------------------
% *** input parameteres:
% *** data - raw EEG  single channel signal
% *** HFOobj.fs - frequency sampling rate
% *** HFOobj.filter - precomputed FIR 
% *** HFOobj.lab_bip - name of the channel
% ***  
% -----OUTPUT--------------------------------------------------------------
% *** output parameteres:
% *** HFOobj.filter: FIR coefficients
% *** HFOobj.channel_name: {Nx1 cell}
% *** HFOobj.fs: sampling frequency
% *** HFOobj.stageIII_flag: flag for the stage III
% *** HFOobj.Datasetup: struct for the stage III
% *** HFOobj.hp: hp for the filter [Hz]
% *** HFOobj.lp: lp for the filter [Hz]
% *** HFOobj.EVENTS: output of STAGE I. check EVENTSlabel frp description
% *** HFOobj.ch: id of current channel under analysis
% *** HFOobj.BLmu: internally set to 0.85
% *** HFOobj.CDFlevel: internally  set to 0.99
% *** HFOobj.BLst_freq: internally  set depending on HFO band
% *** HFOobj.BLborder:  internally  set depending on HFO band
% *** HFOobj.BLmindist: internally  setdepending on HFO band
% *** HFOobj.dur:  internally  setdepending on HFO band
% *** HFOobj.maxIntervalToJoin:  internally  setdepending on HFO band
% *** HFOobj.EVENTSlabel: {1x12 cell}
% *** HFOobj.lf_bound:  internally  setdepending on HFO band
% *** HFOobj.hf_bound:  internally  setdepending on HFO band
% *** HFOobj.Ampl_bound:  internally  setdepending on HFO band
% *** HFOobj.time_thr:  internally  setdepending on HFO band
% *** HFOobj.THRsss: threshold values for each channel [uV]
% *** HFOobj.THR: last computed threshold
% *** HFOobj.indHighEntr: strcut with samples index for the baseline
% *** HFOobj.b_mean: baseline mean for each channel
% *** HFOobj.b_sd: baseline std for each channel
% *** HFOobj.aboveTHR: obsolete
% *** HFOobj.trig_start: obsolete
% *** HFOobj.trig_end: obsolete
% *** HFOobj.EVENTS2: output of STAGE II. check EVENTSlabel frp description
% *** HFOobj.maxcorr: internally set to 0.8000 (stage III)
% *** HFOobj.N_ch: number of channels
% *** HFOobj.EVENTS3: output of STAGE III. check EVENTSlabel frp description
% *** HFOobj.results(HFOobj.N_ch).results struct rearranged from HFOobj.EVENTS3:  with fields:
%             trig_start: start of the event [samples]
%             aboveTHR  : duration of the event [samples]
%             EnergyLF  : Energy in low freq range
%             EnergyR   : Energy in ripple range
%             EnergyFR  : Energy in FRs range
%             Amplpp    : uVpp
%             PowerTrough : istantanous power of the through [uV]
%             Ftrough   : istantanous freq of the through [Hz]
%             PowmaxFR  : istantanous power of the Peak [uV]
%             fmax_FR   : istantanous freq of the Peak [Hz]
%             THR       : threshold of the current channel [uV]
 
% *** ---------------------------------------------------------------------
% *** for Matlab R14a
% *** version 1.0 (Feb 2014)
% *** (c)  Tommaso Fedele 
% *** email: tommaso.fedele@usz.ch


% =========================================================================


function [HFOobj, results] = UMC_IIIstage_st3(Signal, HFOobj)

  
    % STAGEIII
%     HFOobj.stageIII_flag = input.stageIII_flag;
    HFOobj.maxcorr  = 0.8;
%     if isfield(input,'Datasetup')
%         HFOobj.Datasetup =input.Datasetup;
%     end
    
    %----------------------------------------------------------------------
    % Signal must be [points x channel]
    [~, chd] = min(size(Signal)); 
    if chd == 1
        Signal = Signal';
    end
    HFOobj.N_ch = size(Signal,2);


    % ------PREPROCESSING-------------------------------------------------- 
    % filtering
    [Signalfilt, HFOobj] = filter_and_param(Signal, HFOobj);
    
%     % envelope
%     env = abs(hilbert(Signalfilt));
% 
%     for ch = 1:20,HFOobj.N_ch
%         ch
%         tic
% 
%     % ------STAGE I-------------------------------------------------------- 
%     [HFOobj] = baselinethreshold(Signal(:,ch), Signalfilt(:,ch), env(:,ch), HFOobj);
%     HFOobj   = findEOI(env(:,ch),HFOobj,Signalfilt(:,ch));
%     stageItoc = toc
%     tic
%     % ------STAGE II------------------------------------------------------- 
%     HFOobj   = events_validation(Signal(:,ch), Signalfilt(:,ch), HFOobj, ch);
%     stageIItoc = toc
%     end
%     
    if size(HFOobj.EVENTS,1)>1
        if HFOobj.stageIII_flag == 1
            disp stageIII
            tic
        % ------STAGE III------------------------------------------------------ 
            HFOobj  = multichannel_validation(Signalfilt, HFOobj);
            stageIIItoc = toc
        end
    end
    [HFOobj, results] = transfer2results(HFOobj);

end

% ========================================================================= 
% =========================================================================
% =========================================================================
% =========================================================================
function [signalfilt, HFOobj] = filter_and_param(Signal, HFOobj)

    % Filter Signal in the range [hp lp]
    
    
    
    
    switch HFOobj.hp
        case 80  % Ripples

            B = HFOobj.filter.Rb;
            A = HFOobj.filter.Ra;
            HFOobj.lf_bound = 60;
            HFOobj.hf_bound = 250;
            HFOobj.Ampl_bound = 500;
            HFOobj.time_thr = 0.02*HFOobj.fs;
            
        case 250    % Fast Ripples

            B = HFOobj.filter.FRb;
            A = HFOobj.filter.FRa;
            HFOobj.lf_bound = 200;
            HFOobj.hf_bound = 450;
            HFOobj.Ampl_bound = 100;
            HFOobj.time_thr = 0.01*HFOobj.fs;
            
        otherwise

           disp 'WHATrUfiltering????'

    end

    signalfilt=filtfilt(B,A,Signal); %zero-phase filtering
    
     
    
end
 
function HFOobj  = multichannel_validation(Signalfilt, HFOobj);
    
    
    
    HFOobj.EVENTS3         = HFOobj.EVENTS2;
    
    Nev = size(HFOobj.EVENTS3,1);
    
    % building the strucutre suitable for STAGE III
    
    if Nev>0

        % if there are events
        for ev = 1:Nev

            ch = HFOobj.EVENTS3(ev,1);
            start = HFOobj.EVENTS3(ev,2);
            dur = HFOobj.EVENTS3(ev,3);


           elong = 50;
           interval = (start  - elong) : (start +dur + elong);

           %check borders
           interval(interval<1) = [];
           interval(interval>length(Signalfilt)) = [];

           % correlation across more than 4 channels

           % the artitfact is psread over at least 4 channels
           corr_flag = 1
           ch
           if corr_flag 
               [R,p] = corrcoef(Signalfilt(interval,:));
               Rcol = R(ch,HFOobj.Datasetup(ch).Dist_ord);
               pcol = p(ch,HFOobj.Datasetup(ch).Dist_ord);
           else
               [R ] = cov(SignalFilt(interval,:));
               Rcol = R(ch,HFOobj.Datasetup(ch).Dist_ord)/max(R(ch,HFOobj.Datasetup(ch).Dist_ord));
           end

           clear LimCorr    

           LimCorr = ones(1,length(Rcol))* HFOobj.maxcorr;
           mindist = find(HFOobj.Datasetup(ch).Dist_val <=1 );
           LimCorr(mindist) = 1;


           if length(find(abs(Rcol)> LimCorr))

                stageIIIresp(ev) = 1; % ciao
            else
                stageIIIresp(ev) = 0; % stay

            end
        end
        
        HFOobj.EVENTS3(find(stageIIIresp),:) = []; 
        
    end
      
   
end

function [HFOobj, results] = transfer2results(HFOobj);


    if length(HFOobj.EVENTS) == 0
        HFOobj.EVENTS2 = [];
        HFOobj.EVENTS3 = [];
         for ch = 1:HFOobj.N_ch
                    results(ch).results(1).trig_start  = -1;
                    results(ch).results(1).aboveTHR    = -1;
                    results(ch).results(1).EnergyR     = -1;
                    results(ch).results(1).EnergyFR    = -1;
                    results(ch).results(1).Amplpp      = -1;
                    results(ch).results(1).PowerTrough = -1;
                    results(ch).results(1).Ftrough     = -1;
                    results(ch).results(1).PowmaxFR    = -1;
                    results(ch).results(1).fmax_FR     = -1;
                    results(ch).results(1).THR         = -1;  
         end
    end
        
    if length(HFOobj.EVENTS) & isfield(HFOobj,'EVENTS2')
        for ch = 1:HFOobj.N_ch

            if HFOobj.stageIII_flag == 0          
                ecco = find(HFOobj.EVENTS2(:,1) == ch);  
                HFOobj.EVENTS3 = HFOobj.EVENTS2;
            else  
                ecco = find(HFOobj.EVENTS3(:,1)  == ch);
            end

            if length(ecco) > 0

                for evv = 1:length(ecco)

                    results(ch).results(evv).trig_start  = HFOobj.EVENTS3(ecco(evv),2);
                    results(ch).results(evv).aboveTHR    = HFOobj.EVENTS3(ecco(evv),3);
                    results(ch).results(evv).EnergyLF    = HFOobj.EVENTS3(ecco(evv),4);
                    results(ch).results(evv).EnergyR     = HFOobj.EVENTS3(ecco(evv),5);
                    results(ch).results(evv).EnergyFR    = HFOobj.EVENTS3(ecco(evv),6);
                    results(ch).results(evv).Amplpp      = HFOobj.EVENTS3(ecco(evv),7);
                    results(ch).results(evv).PowerTrough = HFOobj.EVENTS3(ecco(evv),8);
                    results(ch).results(evv).Ftrough     = HFOobj.EVENTS3(ecco(evv),9);
                    results(ch).results(evv).PowmaxFR    = HFOobj.EVENTS3(ecco(evv),10);
                    results(ch).results(evv).fmax_FR     = HFOobj.EVENTS3(ecco(evv),11);
                    results(ch).results(evv).THR         = HFOobj.EVENTS3(ecco(evv),12);

                end

            else

                    results(ch).results(1).trig_start  = -1;
                    results(ch).results(1).aboveTHR    = -1;
                    results(ch).results(1).EnergyR     = -1;
                    results(ch).results(1).EnergyFR    = -1;
                    results(ch).results(1).Amplpp      = -1;
                    results(ch).results(1).PowerTrough = -1;
                    results(ch).results(1).Ftrough     = -1;
                    results(ch).results(1).PowmaxFR    = -1;
                    results(ch).results(1).fmax_FR     = -1;
                    results(ch).results(1).THR         = -1;      

            end
        end
    else
        for ch = 1:HFOobj.N_ch
                    results(ch).results(1).trig_start  = -1;
                    results(ch).results(1).aboveTHR    = -1;
                    results(ch).results(1).EnergyR     = -1;
                    results(ch).results(1).EnergyFR    = -1;
                    results(ch).results(1).Amplpp      = -1;
                    results(ch).results(1).PowerTrough = -1;
                    results(ch).results(1).Ftrough     = -1;
                    results(ch).results(1).PowmaxFR    = -1;
                    results(ch).results(1).fmax_FR     = -1;
                    results(ch).results(1).THR         = -1;  
        end
    end
end
