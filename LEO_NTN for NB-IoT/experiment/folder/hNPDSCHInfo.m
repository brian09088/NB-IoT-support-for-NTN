% NB-IoT NPDSCH Information
%
% The class provides NPDSCH related information: the number of
% subframes for a NPDSCH (NSF), the number of repetitions for NPDSCH
% transmission (NRep), the transport block size (TBS) and subframe
% pattern for NPDSCH transmission, based on a set of chosen parameters:
% NPDSCHDataType, ISF, IRep, IMCS, SchedulingInfoSIB1, Modulation and
% OperationMode.
%
% NPDSCHDataType indicates whether the NPDSCH transmission carries
% SystemInformationBlockType1-NB (SIB1-NB) or not, and whether the NPDSCH
% carries the broadcast control channel (BCCH) or not. The allowed values
% are 'SIB1NB', 'BCCHNotSIB1NB' and 'NotBCCH'. Note that SIB1-NB belongs to
% BCCH. When NPDSCH carries SIB1-NB, use ISF and SchedulingInfoSIB1 to
% control NSF, NRep and TBS; When NPDSCH does not carry SIB1-NB, use ISF,
% IRep and IMCS to control NSF, NRep and TBS.

%   Copyright 2017-2023 The MathWorks, Inc.

classdef hNPDSCHInfo
    
    properties (Dependent)
        NPDSCHDataType;     % Character vector to describe the type of the data carried by NPDSCH
        ISF;                % Index of the number of subframes for a NPDSCH (0,1,...,6,7)
        IRep;               % Index of the number of repetitions for NPDSCH transmission (0,1,...,14,15)
        IMCS;               % Index of the modulation and coding scheme (0,1,...,11,13)
        SchedulingInfoSIB1; % Value of schedulingInfoSIB1 (0,1,...,10,11)
        Modulation;         % Modulation used by the NPDSCH (QPSK or 16QAM);
        OperationMode;      % Operation mode ('Standalone','Guardband','Inband-SamePCI','Inband-DifferentPCI')
    end

    properties (Access=private)
        pNPDSCHDataType = 'BCCHNotSIB1NB';
        pISF = 1;
        pIRep = 1;
        pIMCS = 11;
        pSchedulingInfoSIB1 = 0;
        pModulation = 'QPSK';
        pOperationMode = 'Inband-DifferentPCI';
        Inband = true;
    end
    
    properties(Dependent,SetAccess = private) 
        BCCH;   % NPDSCH transmission carrying BCCH or not (true or false)
        SIB1NB; % NPDSCH transmission carrying SystemInformationBlockType1-NB or not (true or false)
        ITBS;   % Index of transport block size
        TBS;    % Transport block size
        NSF;    % Number of subframes for a NPDSCH 
        NRep;   % Number of repetitions for NPDSCH transmission
    end
    
    properties (Constant,Access=public)
        NSFTable = getNSFTable(); % TS 36.213 Table 16.4.1.3-1: Number of subframes (NSF) for NPDSCH
        NRepTable = getNRepTable(); % TS 36.213 Table 16.4.1.3-2: Number of repetitions (NRep) for NPDSCH
        TBSTable = getTBSTable(); % TS 36.213 Table 16.4.1.5.1-1: Transport block size (TBS) table
        NRepTableSIB1 = getSchedulingInfoNReptable(); % TS 36.213 Table 16.4.1.3-3: Number of repetitions for NPDSCH carrying SystemInformationBlockType1-NB
        TBSTableSIB1 = getSchedulingInfoTBStable(); % TS 36.213 Table 16.4.1.5.2-1: Transport block size (TBS) table for NPDSCH carrying SystemInformationBlockType1-NB
    end
    
    methods
        
        function param = get.SIB1NB(obj) 
            if strcmp(obj.NPDSCHDataType,'SIB1NB')
                param = true;
            elseif strcmp(obj.NPDSCHDataType,'BCCHNotSIB1NB')
                param = false;
            elseif strcmp(obj.NPDSCHDataType,'NotBCCH')
                param = false;
            end
        end
        
        function param = get.BCCH(obj) 
            if strcmp(obj.NPDSCHDataType,'SIB1NB')
                param = true;
            elseif strcmp(obj.NPDSCHDataType,'BCCHNotSIB1NB')
                param = true;
            elseif strcmp(obj.NPDSCHDataType,'NotBCCH')
                param = false;
            end
        end

        function param = get.ITBS(obj) 
            if obj.SIB1NB == true
                param = obj.SchedulingInfoSIB1;  
            else
                if strcmpi(obj.Modulation,'16QAM')
                    if obj.Inband
                        param = obj.IMCS+11;
                    else
                        param = obj.IMCS+14;
                    end
                else
                    param = obj.IMCS;
                end
            end
        end
        
        function t = get.TBS(obj)
            if obj.SIB1NB
                m = obj.TBSTableSIB1.ITBS==obj.ITBS;
                t = obj.TBSTableSIB1.TBS(m);
            else
                m = (obj.TBSTable.ITBS==obj.ITBS) & (obj.TBSTable.ISF==obj.ISF);
                if any(m)
                    t = obj.TBSTable.TBS(m);
                else
                    t = [];
                end
            end
        end
        
        function n = get.NSF(obj)
            m = obj.NSFTable.ISF==obj.ISF;
            n = obj.NSFTable.NSF(m);
        end
        
        function n = get.NRep(obj)
            if obj.SIB1NB
                m = obj.NRepTableSIB1.SchedulingInfoSIB1==obj.SchedulingInfoSIB1;
                n = obj.NRepTableSIB1.NRep(m);
            else
                m = obj.NRepTable.IRep==obj.IRep;
                n = obj.NRepTable.NRep(m);
                if strcmpi(obj.Modulation,'16QAM')
                    n = 1;
                end
            end
        end
        
    end
        

    methods 
        
        function layout = displaySubframePattern(obj) 
            % Display the subframe pattern for NPDSCH transmission
            txNsfVec = 0:obj.NSF*obj.NRep-1;
            if obj.BCCH == true
                layout = mod(txNsfVec,obj.NSF);
            else
                layout = mod(floor(txNsfVec/min(obj.NRep,4)),obj.NSF);
            end
            figure;
            image(layout*6);
            set(gca,'YTickLabel',[]);  % Remove Y-label
            xlabel('ms (1 ms per subframe)');
            if obj.BCCH == true
                title('NPDSCH Subframe Repetition Pattern (Carrying BCCH)');
            elseif obj.BCCH == false
                title('NPDSCH Subframe Repetition Pattern (Not carrying BCCH)');
            end
            
        end
        
        % Setter associated with the NPDSCHDataType parameter
        function obj = set.NPDSCHDataType(obj,param) 
            obj.pNPDSCHDataType = param;
        end
        % Getter associated with the NPDSCHDataType parameter
        function param = get.NPDSCHDataType(obj)     
            param = obj.pNPDSCHDataType;      
        end
         
        % Setter associated with the ISF parameter
        function obj = set.ISF(obj,param) 
            m = obj.NSFTable.ISF==param;
            if ~any(m)
                error(['ISF is not one of the set {' num2str(obj.NSFTable.ISF.') '}.']);
            end
            obj.pISF = param;
        end
        % Getter associated with the ISF parameter
        function param = get.ISF(obj)     
            param = obj.pISF;      
        end
        
        % Setter associated with the IRep parameter
        function obj = set.IRep(obj,param)
            if obj.SIB1NB == true
                warning('For NPDSCH carrying SIB1-NB, set SchedulingInfoSIB1 to control number of repetitions (NRep)');
            end
            if strcmpi(obj.Modulation,'16QAM') && ~obj.SIB1NB
                warning('For NPDSCH with 16QAM not carrying SIB1-NB, the number of repetitions (NRep) is always 1');
            end
            m = obj.NRepTable.IRep==param;
            if ~any(m)
                error(['IRep is not one of the set {' num2str(obj.NRepTable.IRep.') '}.']);
            end
            obj.pIRep = param;
        end
        % Getter associated with the IRep parameter
        function param = get.IRep(obj)     
            param = obj.pIRep;      
        end
        
        % Setter associated with the IMCS parameter
        function obj = set.IMCS(obj,param) 
            if obj.SIB1NB == true
                warning('For NPDSCH carrying SIB1-NB, set SchedulingInfoSIB1 to control the transport block size (ITBS or TBS)');
            else
                if ~strcmpi(obj.Modulation,'16QAM') && obj.Inband && (param<0 || param>10) % 0 <= ITBS <= 10, ITBS = IMCS
                    error(['For NPDSCH with QPSK not carrying SIB1-NB and the specified operation mode, IMCS is not one of the set {' num2str(0:10) '}.']);
                elseif ~strcmpi(obj.Modulation,'16QAM') && ~obj.Inband && (param<0 || param>13) % 0 <= ITBS <= 13, ITBS = IMCS
                    error(['For NPDSCH with QPSK not carrying SIB1-NB and the specified operation mode, IMCS is not one of the set {' num2str(0:13) '}.']);
                elseif strcmpi(obj.Modulation,'16QAM') && obj.Inband && (param<0 || param>6) % 11 <= ITBS <= 17, ITBS = IMCS+11
                    error(['For NPDSCH with 16QAM not carrying SIB1-NB and the specified operation mode, IMCS is not one of the set {' num2str(0:6) '}.']);
                elseif strcmpi(obj.Modulation,'16QAM') && ~obj.Inband && (param<0 || param>7) % 14 <= ITBS <= 21, ITBS = IMCS+14
                    error(['For NPDSCH with 16QAM not carrying SIB1-NB and the specified operation mode, IMCS is not one of the set {' num2str(0:7) '}.']);
                end
            end
            obj.pIMCS = param;
        end
        % Getter associated with the IMCS parameter
        function param = get.IMCS(obj)     
            param = obj.pIMCS;      
        end
        
        % Setter associated with the SchedulingInfoSIB1 parameter
        function obj = set.SchedulingInfoSIB1(obj,param) 
            if obj.SIB1NB == false
                warning('For NPDSCH not carrying SIB1-NB, set IMCS to control the transport block size and IRep to control the number of repetitions.');
            end
            m = obj.NRepTableSIB1.SchedulingInfoSIB1==param;
            if ~any(m)
                error(['SchedulingInfoSIB1 is not one of the set {' num2str(obj.NRepTableSIB1.SchedulingInfoSIB1.') '}.']);
            end 
            obj.pSchedulingInfoSIB1 = param;
        end
        % Getter associated with the SchedulingInfoSIB1 parameter
        function param = get.SchedulingInfoSIB1(obj)     
            param = obj.pSchedulingInfoSIB1;
        end

        % Setter associated with the Modulation parameter
        function obj = set.Modulation(obj,param)
            if ~any(strcmpi(param,{'QPSK','16QAM'}))
                error('Modulation is not one of the set {QPSK, 16QAM}');
            end
            obj.pModulation = param;
        end
        % Getter associated with the Modulation parameter
        function param = get.Modulation(obj)
            param = obj.pModulation;
        end

        % Setter associated with the OperationMode parameter
        function obj = set.OperationMode(obj,param)
            if ~any(strcmpi(param,{'Standalone','Guardband','Inband-SamePCI','Inband-DifferentPCI'}))
                error('OperationMode is not one of the set {Standalone, Guardband, Inband-SamePCI, Inband-DifferentPCI}');
            end
            obj.pOperationMode = param;
            if any(strcmpi(param,{'Inband-SamePCI','Inband-DifferentPCI'}))
                obj.Inband = true;
            else
                obj.Inband = false;
            end
        end
        % Getter associated with the OperationMode parameter
        function param = get.OperationMode(obj)
            param = obj.pOperationMode;
        end

    end
    
end

function tab = getNSFTable()
    ISF=(0:7).'; 
    NSF=[1:6 8 10].';
    tab = table(ISF,NSF);
    tab.Properties.Description = 'TS 36.213 Table 16.4.1.3-1: Number of subframes (NSF) for NPDSCH';
end

function tab = getNRepTable()
    IRep=(0:15).'; 
    NRep=[1 2 4 8 16 32 64 128 192 256 384 512 768 1024 1536 2048].';
    tab = table(IRep,NRep);
    tab.Properties.Description = 'TS 36.213 Table 16.4.1.3-2: Number of repetitions (NRep) for NPDSCH';
end

function tab = getTBSTable()
    ITBS=[zeros(1,8)   ones(1,8)   2*ones(1,8) 3*ones(1,8) 4*ones(1,8)  5*ones(1,8) ...
          6*ones(1,8)  7*ones(1,8) 8*ones(1,8) 9*ones(1,8) 10*ones(1,8) 11*ones(1,8) ...
          12*ones(1,8) 13*ones(1,8) 14*ones(1,8) 15*ones(1,8) 16*ones(1,8) 17*ones(1,8)...
          18*ones(1,8) 19*ones(1,8) 20*ones(1,8) 21*ones(1,8)].';
    ISF=repmat((0:7).',[22 1]);
    TBS=[16	 32	 56	  88  120  152  208  256 ...
         24	 56	 88	 144  176  208  256  344 ...
         32	 72	144	 176  208  256  328  424 ...
         40	104	176	 208  256  328  440  568 ...
         56	120	208	 256  328  408  552  680 ...
         72	144	224	 328  424  504  680  872 ...
         88	176	256	 392  504  600  808 1032 ...
        104	224	328	 472  584  680  968 1224 ...
        120	256	392	 536  680  808 1096 1352 ... 
        136	296	456	 616  776  936 1256 1544 ... 
        144	328	504	 680  872 1032 1384 1736 ... 
        176	376	584	 776 1000 1192 1608 2024 ...  
        208	440	680	 904 1128 1352 1800 2280 ... 
        224 488 744 1032 1256 1544 2024 2536 ...
        256 552 840 1128 1416 1736 2280 2856 ...
        280 600 904 1224 1544 1800 2472 3112 ...
        296 632 968 1288 1608 1928 2600 3240 ...
        336 696 1064 1416 1800 2152 2856 3624 ...
        376 776 1160 1544 1992 2344 3112 4008 ...
        408 840 1288 1736 2152 2600 3496 4264 ...
        440 904 1384 1864 2344 2792 3752 4584 ...
        488 1000 1480 1992 2472 2984 4008 4968].';
    tab = table(ITBS,ISF,TBS);
    tab.Properties.Description = 'TS 36.213 Table 16.4.1.5.1-1: Transport block size (TBS) table';
end

function tab = getSchedulingInfoNReptable()
    SchedulingInfoSIB1=(0:11).'; NRep=[4 8 16 4 8 16 4 8 16 4 8 16].';
    tab = table(SchedulingInfoSIB1,NRep);
    tab.Properties.Description = 'TS 36.213 Table 16.4.1.3-3: Number of repetitions for NPDSCH carrying SystemInformationBlockType1-NB';
end

function tab = getSchedulingInfoTBStable()
    ITBS=(0:11).'; 
    TBS=[208 208 208 328 328 328 440 440 440 680 680 680].';
    tab = table(ITBS,TBS);
    tab.Properties.Description = 'TS 36.213 Table 16.4.1.5.2-1: Transport block size (TBS) table for NPDSCH carrying SystemInformationBlockType1-NB';
end
