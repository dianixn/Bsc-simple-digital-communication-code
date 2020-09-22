% OFDM_Base
% Not include training sequence and signal field
% Include Pilot/Data field/CP
% Data can vary even without channel effect and noise
% AWGN Noise is applied to the radiation signal
% MATLAB Build-in OFDM function perform not good because of the absence of
% preamble

SNR = 0;
MSE = 0;
M = 4;
k = log2(M);

%%% Setup OFDM Modulator

mod = comm.OFDMModulator('NumGuardBandCarriers',[4;3],...
'PilotInputPort',true, ...
'PilotCarrierIndices',[12 11; 26 27; 40 39; 54 55], ...%
'NumSymbols',2, ...
'InsertDCNull',true);

%%%

%%% Data generation

modDim = info(mod); % The structure that consists of the dimnesion property in mod
%dataIn = complex(randn(modDim.DataInputSize),randn(modDim.DataInputSize));
dataIn = randi([0 1], modDim.DataInputSize(1), k);
dataSym = bi2de(dataIn);
txSig = pskmod(dataSym,M);
txSig = [txSig, txSig];

pilotIn = complex(rand(modDim.PilotInputSize),rand(modDim.PilotInputSize));
%%%

%%% OFDM Modulate

modData = step(mod,txSig,pilotIn); % Data/CP/Pilot no preamble

%%%

%%% Demodulation setup

demod = comm.OFDMDemodulator(mod);

%%%

%%% Demodulate OFDM Pilot/Data

[dataOut, pilotOut] = step(demod,awgn(modData, SNR));

RxSig1 = pskdemod(dataOut(:,1),M);
RxSig2 = pskdemod(dataOut(:,2),M);
%%% Evaluation

for i = 1 : size(dataIn, 1)
    
    MSE_Individual = norm(dataOut(i, :) - dataIn(i, :), 2);
    MSE = MSE + MSE_Individual;
end

MSE = MSE / size(dataIn, 1);

% BER calculation
numErrs = sum(sum(round(dataSym) ~= round(RxSig1)));
BER = numErrs / length(dataSym);

%%% DataisSameorNot = (max(abs([dataIn(:) - dataOut(:); pilotIn(:) - pilotOut(:)])) < 1e-7)
%isSame = (max(abs([dataIn(:) - dataOut(:); pilotIn(:) - pilotOut(:)])) < 1e-10)