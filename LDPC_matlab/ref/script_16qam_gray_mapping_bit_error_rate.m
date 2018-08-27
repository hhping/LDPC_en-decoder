%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com 
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna Pillai
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 05 June 2008
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bit Error Rate for 16-QAM modulation using Gray modulation mapping

clear
N = 10^5; % number of symbols
M = 16;   % constellation size
k = log2(M); % bits per symbol

% defining the real and imaginary PAM constellation
% for 16-QAM
alphaRe = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
alphaIm = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
k_16QAM = 1/sqrt(10);

Eb_N0_dB  = [0:15]; % multiple Es/N0 values
Es_N0_dB  = Eb_N0_dB + 10*log10(k);

% Mapping for binary <--> Gray code conversion
ref = [0:k-1];
map = bitxor(ref,floor(ref/2));
[tt ind] = sort(map);                                

for ii = 1:length(Eb_N0_dB)
    
    % symbol generation
    % ------------------
    ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
    ipBitReshape = reshape(ipBit,k,N).';
    bin2DecMatrix = ones(N,1)*(2.^[(k/2-1):-1:0]) ; % conversion from binary to decimal
    % real
    ipBitRe =  ipBitReshape(:,[1:k/2]);
    ipDecRe = sum(ipBitRe.*bin2DecMatrix,2);
    ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
    % imaginary
    ipBitIm =  ipBitReshape(:,[k/2+1:k]);
    ipDecIm = sum(ipBitIm.*bin2DecMatrix,2);
    ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2)); 
    % mapping the Gray coded symbols into constellation
    modRe = alphaRe(ipGrayDecRe+1);
    modIm = alphaIm(ipGrayDecIm+1);
    % complex constellation
    mod = modRe + j*modIm;
    s = k_16QAM*mod; % normalization of transmit power to one 
    
    % noise
    % -----
    n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance 
    
    y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    % ------------
    y_re = real(y)/k_16QAM; % real part
    y_im = imag(y)/k_16QAM; % imaginary part

    % rounding to the nearest alphabet
    ipHatRe = 2*floor(y_re/2)+1;
    ipHatRe(find(ipHatRe>max(alphaRe))) = max(alphaRe);
    ipHatRe(find(ipHatRe<min(alphaRe))) = min(alphaRe);
    ipHatIm = 2*floor(y_im/2)+1;
    ipHatIm(find(ipHatIm>max(alphaIm))) = max(alphaIm);
    ipHatIm(find(ipHatIm<min(alphaIm))) = min(alphaIm);

    % Constellation to Decimal conversion
    ipDecHatRe = ind(floor((ipHatRe+4)/2+1))-1; % LUT based
    ipDecHatIm = ind(floor((ipHatIm+4)/2+1))-1; % LUT based

    % converting to binary string
    ipBinHatRe = dec2bin(ipDecHatRe,k/2);
    ipBinHatIm = dec2bin(ipDecHatIm,k/2);

    % converting binary string to number
    ipBinHatRe = ipBinHatRe.';
    ipBinHatRe = ipBinHatRe(1:end).';
    ipBinHatRe = reshape(str2num(ipBinHatRe).',k/2,N).' ;
    
    ipBinHatIm = ipBinHatIm.';
    ipBinHatIm = ipBinHatIm(1:end).';
    ipBinHatIm = reshape(str2num(ipBinHatIm).',k/2,N).' ;

    % counting errors for real and imaginary
    nBitErr(ii) = size(find([ipBitRe- ipBinHatRe]),1) + size(find([ipBitIm - ipBinHatIm]),1) ;

end 
simBer = nBitErr/(N*k);
theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));

close all; figure
semilogy(Eb_N0_dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
axis([0 15 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 16-QAM modulation')








