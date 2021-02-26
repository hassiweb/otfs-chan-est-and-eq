% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef LteChannelCoder < matlab.System
    
    properties
        % Input properties
        TurboDecodingNumIterations;
        ModOrder;
        NumLayers;
        OutputLength;
        LinkDirection;
        RedundancyVersion;
        TBS;
        
        % Init
        L;
        C;
        Cplus;
        Kplus;
        Kminus;
        F;

        intrlvrIndicesForKplus;
        intrlvrIndicesForKminus;
        
        Crc24ACoder;
        Crc24BCoder;
        LteTurboCoder;
        
        % For debug
%         bitSeqA;
%         bitSeqB;
%         bitSeqC;
%         bitSeqD;
%         bitSeqE;
%         bitSeqF;
    end
    
    methods
        %% Constructor
        function this = LteChannelCoder(varargin) % Constructor
            if nargin > 0
                setProperties(this, nargin, varargin{:}, 'TurboDecodingNumIterations', 'ModOrder', 'NumLayers', 'OutputLength', 'LinkDirection', 'RedundancyVersion', 'TBS');

                % This flow is defined in TS36.212 section 5.1.2
                B = this.TBS + 24; % input bit length (TBS + CRC24bits)
                Z = 6144;
                if (B <= Z)
                    L = 0; % CRC(gCRC24B) length
                    C = 1; % number of code blocks
                    Bprime = B; % total bits
                else
                    L = 24; % CRC length
                    C = ceil(B/(Z-L)); % number of code blocks
                    Bprime = B + C*L; % total bits
                end
                % Values of block size K from table 5.1.3-3
                Ktable = [40:8:512 528:16:1024 1056:32:2048 2112:64:6144].';
                % First segmentation size: Kplus
                temp = find(Ktable >= Bprime/C);
                Kplus = Ktable(temp(1), 1); % minimum K in Ktable such that C*K >= Bprime

                if C == 1
                    Kminus = 0;
                    Cplus  = 1;
                    Cminus = 0;
                else
                    % Second segmentation size: Kminus
                    temp2 = find(Ktable < Kplus);
                    Kminus = Ktable(temp2(end), 1); % maximum K in Ktable such that K < Kplus
                    Cminus = floor((C*Kplus-Bprime)/(Kplus-Kminus));
                    Cplus  = C - Cminus;
                end

                % Generate segmented bit sequence
                F = Cplus * Kplus + Cminus * Kminus - Bprime;

                % Store variables
                this.L = L;
                this.C = C;
                this.Cplus = Cplus;
                this.Kplus = Kplus;
                this.Kminus = Kminus;
                this.F = F;

                % Turbo code internal interleaver, defined in TS36.212 section 5.1.3.2.3
                % For Kplus
                [f1, f2] = getf1f2(Kplus); % defined in TS36.212 section 5.1.3.2.3 Table 5.1.3-3
                Idx      = (0:Kplus-1).';
                this.intrlvrIndicesForKplus = mod(f1*Idx + f2*Idx.^2, Kplus) + 1;
                % For Kminus
                if C > 1
                    [f1, f2] = getf1f2(Kminus); % defined in TS36.212 section 5.1.3.2.3 Table 5.1.3-3
                    Idx      = (0:Kminus-1).';
                    this.intrlvrIndicesForKminus = mod(f1*Idx + f2*Idx.^2, Kminus) + 1;
                end

                % CRC 24 bit
                this.Crc24ACoder = Crc24Coder('A');
                this.Crc24BCoder = Crc24Coder('B');
                % Turbo encoder/decoder
                this.LteTurboCoder = TurboCoder( ...
                    'InterleaverIndicesSource', 'Input port', ...
                    'NumIterations', this.TurboDecodingNumIterations);
            end
        end
        
        
        
        %% Encoder
        function bitSeqF = encode(this, in)
            bitSeqA = in;
            bitSeqB = this.Crc24ACoder.generate(bitSeqA);
            bitSeqC = this.processCodeBlockSegmentation(bitSeqB);
            bitSeqD = cell(this.C,1); bitSeqE = cell(this.C,1);
            for r=0:this.C-1
                bitSeqD{r+1} = this.processLteTurboEncoding(bitSeqC{r+1}, r);
                bitSeqE{r+1} = this.processRateMatching(bitSeqD{r+1}, r);
            end
            bitSeqF = this.processCodeBlockConcatenation(bitSeqE);
%             this.bitSeqA=bitSeqA; this.bitSeqB=bitSeqB; this.bitSeqC=bitSeqC; this.bitSeqD=bitSeqD; this.bitSeqE=bitSeqE; this.bitSeqF=bitSeqF;
        end
        
        %% Decoder
        function [bitSeqA, error] = decode(this, in)
            bitSeqF = in;
            bitSeqE = this.processCodeBlockSeparation(bitSeqF);
            bitSeqD = cell(this.C,1); bitSeqC = cell(this.C,1);
            for r=0:this.C-1
                bitSeqD{r+1} = this.processRateRecovering(bitSeqE{r+1}, r);
                bitSeqC{r+1} = this.processLteTurboDecoding(bitSeqD{r+1}, r);
            end
            [bitSeqB, blkErrB] = this.processCodeBlockDesegmentation(bitSeqC);
            [bitSeqA, blkErrA] = this.Crc24ACoder.detect(bitSeqB);
            error = blkErrA | blkErrB;
        end
        
        
        
        %% Internal functions
        % Code block segmentation
        function bitSeqC = processCodeBlockSegmentation(this, bitSeqB)
            bitSeqC = cell(this.C,1);
            for r = 0:this.C-1
                if r < this.Cplus
                    Kr = this.Kplus;
                else
                    Kr = this.Kminus;
                end
                bitSeqC{r+1} = zeros(Kr, 1);
            end
            for k=0:this.F-1
                bitSeqC{1}(k+1) = nan;
            end
            k = this.F;
            s = 0;
            for r = 0:this.C-1
                if r < this.Cplus
                    Kr = this.Kplus;
                else
                    Kr = this.Kminus;
                end
                while k < Kr-this.L
                    bitSeqC{r+1}(k+1)=bitSeqB(s+1);
                    k = k+1;
                    s = s+1;
                end
                if this.C>1
                    bitSeqC{r+1}(:) = this.Crc24BCoder.generate(bitSeqC{r+1}(1:Kr-this.L));
                end
                k=0;
            end
        end
        
        function [bitSeqB, blkErr] = processCodeBlockDesegmentation(this, bitSeqC)
            if this.F > 0 % Remove dummy bits
                bitSeqC{1}(1:this.F) = [];
            end

            blkErr = false;            
            bitSeqB = zeros(length(cell2mat(bitSeqC)),1)*NaN;
            index = 1;
            if this.C > 1
                for r=0:this.C-1
                    % CRC detector
                    [bitSeqBtmp, blkErrTmp] = this.Crc24BCoder.detect(bitSeqC{r+1});
                    bitLength = length(bitSeqBtmp);
                    bitSeqB(index:index+bitLength-1) = bitSeqBtmp;
                    index = index + bitLength;
                    blkErr = blkErr || blkErrTmp;
                end
                bitSeqB(isnan(bitSeqB)) = [];
            else
                bitSeqB = bitSeqC{1};
            end
        end

        % LTE Turbo encoding
        function bitSeqD = processLteTurboEncoding(this,bitSeqC, r)
            % K is the number of bits to encoding
            if r < this.Cplus
                K = this.Kplus;
                intrlvrIndices = this.intrlvrIndicesForKplus;
            else
                K = this.Kminus;
                intrlvrIndices = this.intrlvrIndicesForKminus;
            end

            % If the code block to be encoded is the 0-th code block and the number of filler bits is greater than zero, i.e., F > 0,
            % then the encoder shall set ck, = 0, k = 0,...,(F-1) at its input
            if r == 0 && this.F > 0
                bitSeqC(~isfinite(bitSeqC)) = 0; % Replacement of Dummy Bits
            end

            % Turbo encoding
            bitSeqD = this.LteTurboCoder.encode(bitSeqC(1:K), intrlvrIndices);

            % If the code block to be encoded is the 0-th code block and the number of filler bits is greater than zero, i.e., F > 0,
            % then the encoder shall set d0=<NULL>, k = 0,...,(F-1) and d1=<NULL>, k = 0,..,(F-1) at its output.
            if r == 0 && this.F > 0
                % Bit streams
                d0 = bitSeqD(1:3:end); % systematic
                d1 = bitSeqD(2:3:end); % parity 1st
                d0(1:this.F) = nan;
                d1(1:this.F) = nan;
                bitSeqD(1:3:end) = d0;
                bitSeqD(2:3:end) = d1;
            end
        end
        
        function bitSeqC = processLteTurboDecoding(this, bitSeqD, r)
            % Turbo decoding
            % Replace filter bits with zero %EDIT: should not be zero LLR since
            % we know their value: it should be Inf
            if r == 0 && this.F > 0
                % Bit streams
                d0 = bitSeqD(1:3:end); % systematic
                d1 = bitSeqD(2:3:end); % parity 1st
%                 d0(1:this.F) = -Inf;
%                 d1(1:this.F) = -Inf;
                d0(1:this.F) = Inf;
                d1(1:this.F) = Inf;
                bitSeqD(1:3:end) = d0;
                bitSeqD(2:3:end) = d1;
            end        

            % Calculate interleaver index in the turbo decoder
            if r < this.Cplus
                intrlvrIndices = this.intrlvrIndicesForKplus;
            else
                intrlvrIndices = this.intrlvrIndicesForKminus;
            end

            % Avoid infinity calculation
            bitSeqD(bitSeqD==Inf) = 700;
            bitSeqD(bitSeqD==-Inf) = -700;
            
            bitSeqC = this.LteTurboCoder.decode(-bitSeqD, intrlvrIndices);
%             [~,~,bitSeqC,~] = this.LteTurboCoder.decodeWithUpdatedLlrOutput(-bitSeqD, intrlvrIndices, this.Crc24ACoder);
        end
        
        % Rate matching per coded block, with and without the bit selection.
        function bitSeqE = processRateMatching(this, bitSeqD, r)
            % K is the number of bits to encoding
            if r < this.Cplus
                K = this.Kplus;
            else
                K = this.Kminus;
            end
            D = K + 4; % D is the number of encoded bits per output stream

            if numel(bitSeqD)~=3*D, error('D times 3 must be size of input 1.');end

            % Parameters
            colTcSb = 32;              % number of columns of the sub-block interleaver matrix, which is defined in 3GPP TS36.211 section 5.1.4.1.1 
            rowTcSb = ceil(D/colTcSb); % number of rows of the sub-block interleaver matrix
            Kpi = colTcSb * rowTcSb;   % number of elements of the matrix, number of interleaved bits of each stream
            Nd  = Kpi - D;             % number of dummy bits for padding

            % Bit streams
            d0 = bitSeqD(1:3:end); % systematic
            d1 = bitSeqD(2:3:end); % parity 1st
            d2 = bitSeqD(3:3:end); % parity 2nd

            i0 = (1:D)';
            Index  = this.generateSubblkIntrlvrIndices( i0, colTcSb, rowTcSb, Nd); % for d0 and d1 streams
            Index2 = this.generateSubblkIntrlvrIndices2(i0, colTcSb, rowTcSb, Nd); % for d2 stream

            % Sub-block interleaving - per stream, as defined in section 5.1.4.1.1 of TS36.212
            v0 = this.processSubblkIntrlv(d0,Index);
            v1 = this.processSubblkIntrlv(d1,Index);
        %     v2 = this.processSubblkIntrlv(d2,Index);
            v2 = this.processSubblkIntrlv(d2,Index2);

            % Concat 0, interleave 1, 2 sub-blk streams
            vpre = [v1,v2].';
            v12  = vpre(:);

            wk            = zeros(3*Kpi, 1);
            wk(1:Kpi)     = v0;
            wk(Kpi+1:end) = v12;

            % Bit collection, selection as defined in section 5.1.4.1.2 of TS36.212
            % Calculate Ncb which is the soft buffer size and is used to do rate matching
            Kw      = 3*Kpi;   % circular buffer length
            Nsoft   = 3667200;%250368; % the total number of soft channel bits defined in TS36.304 according to the UE category
            Kc      = 1;
            Kmimo   = 1; % 1 for TM1,2,5,6,7,10; 2 for TM3,4,8,9
            Mlimit  = 8; % constant
            Mdlharq = 8; % the maximum number of DL HARQ processes as defined in section 7 of TS36.213
            Nir     = floor(Nsoft/(Kc*Kmimo*min(Mdlharq, Mlimit))); % the soft buffer size for the transport block by Nir bits
            switch this.LinkDirection
                case 'Downlink'
                    Ncb = min( floor(Nir/this.C), Kw ); % the soft buffer size for the r-th code block by Ncb bits
                case 'Uplink'
                    Ncb = Kw;
            end

            % Calculate E which is the rate matching output sequence length for the r-th coded block
            Nl     = this.NumLayers; % 2 for transmit diversity, otherwise, number of layers a transport block is mapped onto
            Qm     = log2(this.ModOrder); % 2 for QPSK, 4 for 16QAM and 6 for 64QAM
            G      = this.OutputLength; % total number of bits available for the transmission of one transport block
            Gprime = G/(Nl*Qm);
            gamma  = mod(Gprime, this.C);

            if r <= this.C-gamma-1
                E = Nl*Qm*floor(Gprime/this.C);
            else
                E = Nl*Qm*ceil(Gprime/this.C);
            end
            k0 = rowTcSb*(2*ceil(Ncb/(8*rowTcSb))*this.RedundancyVersion+2);

            % Create rate matching output sequence
            k = 0;
            j = 0;
            bitSeqE = false(E,1);
            while k < E
                if ~isnan(wk(mod(k0+j,Ncb)+1))
                    bitSeqE(k+1) = wk(mod(k0+j,Ncb)+1);
                    k = k+1;
                end
                j = j+1;
            end
        end
        
        function bitSeqD = processRateRecovering(this, bitSeqE, r)
            % Parameters
            % Define K of current code block
            if r < this.Cplus
                K = this.Kplus;
            else
                K = this.Kminus;
            end
            D       = K+4;
            colTcSb = 32;
            rowTcSb = ceil(D/colTcSb);
            Kpi     = colTcSb * rowTcSb;
            Nd      = Kpi - D;

            Kw      = 3*Kpi;   % circular buffer length
            Nsoft   = 3667200; %250368; % the total number of soft channel bits defined in TS36.304 according to the UE category
            Kc      = 1;
            Kmimo   = 1; % 1 for TM1,2,5,6,7,10; 2 for TM3,4,8,9
            Mlimit  = 8; % constant
            Mdlharq = 8; % the maximum number of DL HARQ processes as defined in section 7 of TS36.213
            Nir     = floor(Nsoft/(Kc*Kmimo*min(Mdlharq, Mlimit))); % the soft buffer size for the transport block by Nir bits
            switch this.LinkDirection
                case 'Downlink'
                    Ncb = min( floor(Nir/this.C), Kw ); % the soft buffer size for the r-th code block by Ncb bits
                case 'Uplink'
                    Ncb = Kw;
            end
            k0 = rowTcSb*(2*ceil(Ncb/(8*rowTcSb))*this.RedundancyVersion+2); % offset

            % Recreate buffer at RX
            %   Sub-block interleaving - per stream - for NAN location in buffer
            i0 = (1:D)';
            Index  = this.generateSubblkIntrlvrIndices( i0, colTcSb, rowTcSb, Nd);
            Index2 = this.generateSubblkIntrlvrIndices2(i0, colTcSb, rowTcSb, Nd);

            d01prime = zeros(D,1);
            if r == 0 && this.F > 0
                d01prime(1:this.F) = NaN;
            end
            d2prime = zeros(D,1);

            v01prime = this.processSubblkIntrlv(d01prime,Index);
        %     v2prime  = this.processSubblkIntrlv(d2prime, Index);
            v2prime  = this.processSubblkIntrlv(d2prime, Index2);

            %   Concat 0, interleave 1, 2 sub-blk streams
            vpreprime = [v01prime,v2prime].';
            v12prime  = vpreprime(:);

            %   Bit collection
            wk            = zeros(3*Kpi, 1);
            wk(1:Kpi)     = v01prime;
            wk(Kpi+1:end) = v12prime; % has the NANs at the right locations

            %   Fill incoming data at the right location accounting for offset (k0)
            %   and NaNs in wk (using a circular buffer). 
            %   Trailing punctures are already zero'ed in the buffer.
            % Calculate E which is the rate matching output sequence length for the r-th coded block
            Nl     = this.NumLayers; % 2 for transmit diversity, otherwise, number of layers a transport block is mapped onto
            Qm     = log2(this.ModOrder); % 2 for QPSK, 4 for 16QAM and 6 for 64QAM
            G      = this.OutputLength; % total number of bits available for the transmission of one transport block
            Gprime = G/(Nl*Qm);
            gamma  = mod(Gprime, this.C);
            if r <= this.C-gamma-1
                E = Nl*Qm*floor(Gprime/this.C);
            else
                E = Nl*Qm*ceil(Gprime/this.C);
            end
            
            k = 0;
            j = 0;
            while k < E
                if ~isnan(wk(mod(k0+j,Ncb)+1))
                    wk(mod(k0+j,Ncb)+1) = bitSeqE(k+1);
                    k = k+1;
                end
                j = j+1;
            end

            % Inverse bit collection - create bit streams
            v0 = wk(1:Kpi);       % systematic
            v1 = wk(Kpi+1:2:end); % parity 1st
            v2 = wk(Kpi+2:2:end); % parity 2nd    

            d0 = this.processSubblkDeintrlv(v0, Index);
            d1 = this.processSubblkDeintrlv(v1, Index);
            d2 = this.processSubblkDeintrlv(v2, Index2);

            % Concat 1, 2, 3 streams - for turbo decoding
            temp = [d0 d1 d2].';
            bitSeqD = temp(:);            
        end

        
        function v = generateSubblkIntrlvrIndices(~, d, colTcSb, rowTcSb, Nd)
            % Sub-block interleaving - for d0 and d1 streams only
            colPermPat = [0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30,...
                          1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31];
            % For 1 and 2nd streams only
            y          = [NaN*ones(Nd, 1); d];       % null (NaN) filling
            inpMat     = reshape(y, colTcSb, rowTcSb).';
            permInpMat = inpMat(:, colPermPat+1);
            v          = permInpMat(:);
        end
        
        function v = generateSubblkIntrlvrIndices2(~, d, colTcSb, rowTcSb, Nd)
            % Sub-block interleaving - for d2 stream only
            K_pi = colTcSb*rowTcSb;
            colPermPat = [0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30,...
                          1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31];
            pi = mod(colPermPat(floor((0:K_pi-1)/rowTcSb)+1) + colTcSb*(mod(0:K_pi-1, rowTcSb)) + 1,K_pi)';
            % For 3rd stream only
            y      = [NaN*ones(Nd, 1); d];       % null (NaN) filling
            inpMat = reshape(y, colTcSb, rowTcSb).';
            ytemp  = inpMat.';
            y      = ytemp(:);
            v      = y(1+pi);
        end

        function out = processSubblkIntrlv(~, d0,Index)
            out         = zeros(size(Index));
            IndexG      = find(~isnan(Index)==1);
            IndexB      = find(isnan(Index)==1);
            out(IndexG) = d0(Index(IndexG));
            Nd          = numel(IndexB);
            out(IndexB) = nan*ones(Nd,1);
        end
        
        function out = processSubblkDeintrlv(~, v0,Index)
            dataIdx         = find(isnan(Index)==0);
            intrlvData      = v0(dataIdx);
            tmpVec          = (1:numel(dataIdx)).';
            intrlvTmpVec    = tmpVec(Index(dataIdx));
            [~,deintrlvIdx] = sort(intrlvTmpVec);
            out             = intrlvData(deintrlvIdx);
        end

        function bitSeqF = processCodeBlockConcatenation(this, bitSeqE)
            k = 0;
            r = 0;
            bitSeqF = false(this.OutputLength,1);
            while r < this.C
                % 3GPP definition
%                 j = 0;
%                 e_rk = bitSeqE{r+1};
%                 while j < numel(e_rk)
%                     bitSeqF(k+1) = e_rk(j+1);
%                     k = k+1;
%                     j = j+1;
%                 end
                bitSeqF(k+1:k+length(bitSeqE{r+1}))= bitSeqE{r+1};
                k = k + length(bitSeqE{r+1});
                r = r+1;
            end
        end
        
        function bitSeqE = processCodeBlockSeparation(this, bitSeqF)
            % Calculate E which is the rate matching output sequence length for the r-th coded block
            Nl     = this.NumLayers; % 2 for transmit diversity, otherwise, number of layers a transport block is mapped onto
            Qm     = log2(this.ModOrder); % 2 for QPSK, 4 for 16QAM and 6 for 64QAM
            G      = this.OutputLength; % total number of bits available for the transmission of one transport block
            Gprime = G/(Nl*Qm);
            gamma  = mod(Gprime, this.C);
            
            bitSeqE = cell(this.C,1);
            currentIdx = 1;            
            for r=0:this.C-1
                if r <= this.C-gamma-1
                    E = Nl*Qm*floor(Gprime/this.C);
                else
                    E = Nl*Qm*ceil(Gprime/this.C);
                end

                bitSeqE{r+1} = bitSeqF(currentIdx:currentIdx+E-1);
                currentIdx = currentIdx + E;            
            end
        end
    end
end