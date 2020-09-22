%default as Gray coding

M = 4;                 % Modulation order
k = log2(M);            % Bits per symbol
EbNoVec = (-5 : 1 : 10)';      % Eb/No values (dB)
numSymPerFrame = 1000;   % Number of QAM symbols per frame

berEst = zeros(size(EbNoVec)); %Initialize the results vector
SNR_dB = zeros(size(EbNoVec));

for n = 1:length(EbNoVec)
    
    % Convert Eb/No to SNR
    SNR_dB(n) = EbNoVec(n) + 10 * log10(k);
    
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;
    
    while numErrs < 200 && numBits < 1e7
        % Generate binary data and convert to symbols
        dataIn = randi([0 1],numSymPerFrame,k);
        %x = randi([0,1],1,N); %Generate random signal
        dataSym = bi2de(dataIn);
        
        % QAM modulate using 'Gray' symbol mapping
        txSig = pskmod(dataSym,M);
        
        % Pass through AWGN channel
        rxSig = awgn(txSig,SNR_dB(n),'measured');
        
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

[berTheory, serTheory] = berawgn(EbNoVec, 'psk', M, 'nondiff'); %Theory value from the build-in funcation

ber_QPSK = qfunc(sqrt( 10.^(SNR_dB/10))); %Theory value from the formula

figure
semilogy(SNR_dB,berEst,'yo')
hold on
semilogy(SNR_dB,berTheory, 'r*')
hold on
semilogy(SNR_dB,ber_QPSK, 'b-')
grid on
legend('Estimated BER','Theoretical BER from berawgn','Theoretical BER')
xlabel('SNR (dB)')
ylabel('Bit Error Rate')

%figure
%semilogy(EbNoVec,berEst,':go')
%hold on
%semilogy(EbNoVec,berTheory, 'y*')
%hold on
%semilogy(EbNoVec,ber_QPSK, 'b--o')
%grid on
%legend('Estimated BER','Theoretical BER','Theoretical BER from formula')
%xlabel('Eb/No in dB')
%ylabel('Bit Error Rate')