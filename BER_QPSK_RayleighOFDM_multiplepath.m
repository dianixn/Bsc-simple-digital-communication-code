% Must use pathdalay < 1e-5, where cof = [1;1]
% default as Gray coding
% comm.AWGNChannel, comm.RicianChannel, comm.MIMOChannel,
% doppler(doppler(specType, varargin))

M = 4;
k = log2(M);            % Bits per symbol
EbNoVec = (-5 : 1 : 5)';      % Eb/No values (dB)
bitRate = 1e3;    % Data rate is 50 kb/s

berEst = zeros(size(EbNoVec)); %Initialize the results vector

%rayleigh_channel = comm.RayleighChannel('PathDelays',[0 2e-5],'AveragePathGains',[0 -9]);

%   chan = comm.RayleighChannel(...
%       'SampleRate',          10e3,...
%       'PathDelays',          [0 1.5e-4],...
%       'AveragePathGains',    [-2 -3],...
%       'NormalizePathGains',  true,...
%       'MaximumDopplerShift', 30,...
%       'DopplerSpectrum',     {doppler('Gaussian',0.6), doppler('Flat')},...
%       'RandomStream',        'mt19937ar with seed',...
%       'Seed',                22,...
%       'PathGainsOutputPort', true);

OFDM_baseband = comm.OFDMModulator('FFTLength',64,...
'NumGuardBandCarriers',[4;3],...
'CyclicPrefixLength',16, ...
'PilotInputPort',true, ...
'PilotCarrierIndices',[12; 26; 40; 54], ...
'NumSymbols',1, ...
'Windowing',true, ...
'InsertDCNull',true);

modDim = info(OFDM_baseband);
numSymPerFrame = modDim;
pilotIn = complex(rand(modDim.PilotInputSize),rand(modDim.PilotInputSize));

rayleigh_channel = comm.RayleighChannel(...
       'SampleRate',          bitRate,...
       'PathDelays',          [0 1e-5],...
       'AveragePathGains',    [0 -2],...
       'MaximumDopplerShift', 1,...
       'DopplerSpectrum',     {doppler('Gaussian',0.6), doppler('Flat')},...
       'NormalizePathGains',  true,...
       'PathGainsOutputPort', true);

OFDM_Demod = comm.OFDMDemodulator(OFDM_baseband);

for n = 1:length(EbNoVec)
    
    % Convert Eb/No to SNR
    SNR_dB = EbNoVec(n) + 10 * log10(k);
    
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;
    
    while numBits < 1e5
        % Generate binary data and convert to symbols
        dataIn = randi([0 1], modDim.DataInputSize(1), k);
        %x = randi([0,1],1,N); %Generate random signal
        dataSym = bi2de(dataIn);
        
        % QAM modulate using 'Gray' symbol mapping
        txSig = pskmod(dataSym,M);
        
        %qpskMod = comm.QPSKModulator( ...
        %'BitInput',    true, ...
        %'PhaseOffset', pi/4);
        
        %%% Data generation and OFDM Modulate
        
        pilotIn = complex(rand(modDim.PilotInputSize),rand(modDim.PilotInputSize));
        modData = step(OFDM_baseband,txSig,pilotIn);
        
        % Fading channel applied, Rayleigh channel
        [fadedSig, path_gain] = rayleigh_channel(modData); % Apply the channel effects
        %impulse(rayleigh_channel)
        [fadedSig1, h1] = step(rayleigh_channel,modData);
        %
        chaninfo = info(rayleigh_channel);
        coeff = chaninfo.ChannelFilterCoefficients;
        
        Np = length(rayleigh_channel.PathDelays);
        fracdelaydata = zeros(size(modData, 1), Np);
        
        for i = 1 : Np
            fracdelaydata(:,i) = filter(coeff(i,:),1,modData);
        end
        
        fading_signal_multipath = path_gain .* fracdelaydata;
        fadedSig2 = sum(fading_signal_multipath, 2);
        
        isequal(fadedSig,fadedSig2)
        
        % Pass through AWGN channel
        rxSig = awgn(fadedSig,SNR_dB,'measured');
        
        % CSI and equailent
        overall_pathgain = sum(path_gain, 2)/ Np;
        rxSig = rxSig ./ overall_pathgain;
        % Demodulate OFDM Pilot/Data
        
        [rx_OFDM, pilotOut] = step(OFDM_Demod,rxSig);
        
        % Demodulate the noisy signal
        rxSym = pskdemod(rx_OFDM, M);
        % Convert received symbols to bits
        dataOut = de2bi(rxSym,k);
        
        % Calculate the number of bit errors
        nErrors = biterr(dataIn,dataOut);
        
        % Increment the error and bit counters
        numErrs = numErrs + nErrors;
        numBits = numBits + modDim.DataInputSize(1) * k;
    end
    
    % Estimate the BER
    berEst(n) = numErrs/numBits;
end

divorder = 2;
[berTheory, serTheory] = berfading(EbNoVec, 'psk', M, divorder); %Theory value from the build-in funcation

%For cases where diversity is used, the Eb/N0 on each diversity branch is EbNo/divorder,
%where divorder is the diversity order (the number of diversity branches) and is a positive integer.

figure
semilogy(EbNoVec,berEst,'r')
%hold on
%semilogy(EbNoVec,berTheory, 'b')
grid on
legend('Estimated BER')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')