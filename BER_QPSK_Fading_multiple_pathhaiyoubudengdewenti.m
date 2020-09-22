%default as Gray coding
% comm.AWGNChannel, comm.RicianChannel, comm.MIMOChannel,
% doppler(doppler(specType, varargin))

M = 4;
k = log2(M);            % Bits per symbol
EbNoVec = (5 : 1 : 20)';      % Eb/No values (dB)
numSymPerFrame = 100;   % Number of QAM symbols per frame
bitRate = 1e3;    % Data rate is 50 kb/s

berEst = zeros(size(EbNoVec)); %Initialize the results vector

%rayleigh_channel = comm.RayleighChannel('PathDelays',[0 2e-5],'AveragePathGains',[0 -9]);

%   chan = comm.RayleighChannel(...
%       'SampleRate',          10e3,...
%       'PathDelays',          [0 1.5e-4],...
%       'AveragePathGains',    [2 3],...
%       'NormalizePathGains',  true,...
%       'MaximumDopplerShift', 30,...
%       'DopplerSpectrum',     {doppler('Gaussian',0.6), doppler('Flat')},...
%       'RandomStream',        'mt19937ar with seed',...
%       'Seed',                22,...
%       'PathGainsOutputPort', true);

rayleigh_channel = comm.RayleighChannel(...
       'SampleRate',          bitRate,...
       'PathDelays',          [0 1.5e-5],...
       'AveragePathGains',    [-2 -3],...
       'NormalizePathGains',  true,...
       'MaximumDopplerShift', 1,...
       'DopplerSpectrum',     {doppler('Gaussian',0.6), doppler('Flat')},...
       'PathGainsOutputPort', true);
   
for n = 1:length(EbNoVec)
    
    % Convert Eb/No to SNR
    SNR_dB = EbNoVec(n) + 10 * log10(k);
    
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;
    
    while numBits < 1e3 %numErrs < 200 && numBits < 1e7
        % Generate binary data and convert to symbols
        dataIn = randi([0 1],numSymPerFrame,k);
        %x = randi([0,1],1,N); %Generate random signal
        dataSym = bi2de(dataIn);
        
        % QAM modulate using 'Gray' symbol mapping
        txSig = pskmod(dataSym,M);
        
        [fadedSig, path_gain] = rayleigh_channel(txSig); % Apply the channel effects
        
        %
        chaninfo = info(rayleigh_channel);
        coeff = chaninfo.ChannelFilterCoefficients;
        Np = length(rayleigh_channel.PathDelays);
        fracdelaydata = zeros(size(txSig,1),Np);
        i =0;
        for i = 1:Np
            fracdelaydata(:,i) = filter(coeff(i,:),1,txSig);
        end
        
        fading_signal_multipath = path_gain .* fracdelaydata;
        fadedSig2 = sum(fading_signal_multipath, 2);
        
        isequal(fadedSig,fadedSig2)
        
        % Pass through AWGN channel
        %rxSig = awgn(fadedSig,SNR_dB,'measured');
        
        % CSI and equailent
        fracdelaydataRx = zeros(size(txSig,1),Np);
        i = 0;
        for i = 1:Np
            fracdelaydataRx(:,i) = filter(1, coeff(i,:), fadedSig);
            rxSig = fracdelaydataRx(:,i)./path_gain(:, i);
        end
        
        isequal(txSig,rxSig)
        
        % Demodulate the noisy signal
        rxSym = pskdemod(rxSig,M);
        % Convert received symbols to bits
        dataOut = de2bi(rxSym,k);
        
        % Calculate the number of bit errors
        nErrors = biterr(dataIn,dataOut);
        
        % Increment the error and bit counters
        numErrs = numErrs + nErrors;
        numBits = numBits + numSymPerFrame*k;
    end
    
    % Estimate the BER
    berEst(n) = numErrs/numBits;
end

divorder = 1;
[berTheory, serTheory] = berfading(EbNoVec, 'psk', M, divorder); %Theory value from the build-in funcation

%For cases where diversity is used, the Eb/N0 on each diversity branch is EbNo/divorder,
%where divorder is the diversity order (the number of diversity branches) and is a positive integer.

figure
semilogy(EbNoVec,berEst,'r')
hold on
semilogy(EbNoVec,berTheory, 'b')
grid on
legend('Estimated BER','Theoretical BER')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')