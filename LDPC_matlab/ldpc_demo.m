% All rights reserved by Chenzy
% The file simulate CMMB-Stimi LDPC simulation model
% Author      : Chen Zhengyi
% Email       : czhengyi@126.com
% Version     : 1.0
% Data        : 25 Dec 2009
%
clc;
clear all;
%modulation = ['bpsk' 'qpsk' '16qm'];
constellation = [ 1 2 4];

% Parameter Definition here
coderate = 0.5;            % 0.75
mode = 'qpsk';             % 'bpsk' 'qpsk' '16qm'
% snr  =    ;              % Signal-Noise ratio
k = constellation(2);      % 1 2 4
N = 1;                     % 500 frame 10^6

%Eb_N0_dB = [1.2:0.1:2]; % multiple Es/N0 values
Eb_N0_dB = [1.2];
Es_N0_dB = Eb_N0_dB + 10*log10(k);

for ii = 1:length(Eb_N0_dB)
    % symbol generation
    % -----------------
    inbit = randint(1,N*9216*coderate);
    if( coderate == 0.5)
        load data\G.mat
        load data\col_order.mat
    elseif( coderate == 0.75)
        load data\G34.mat
        load data\col_order34.mat
    end
    for nframe = 1 : N
        symbol = inbit((nframe-1)*9216*coderate+1:nframe*9216*coderate);
        pc = mod((G*symbol'),2);
        msg(1:9216*(1-coderate)) = pc';
        msg(9216*(1-coderate)+1:9216) = symbol;
        for i=1:9216
            code(col_order(i)+1)=msg(i);
        end
        tx(1,(nframe-1)*9216+1:nframe*9216)=code;
    end
    % clear inbit;
    clear msg;
    clear symbol;
    clear pc;
    clear code;
    clear G;
    % Complex constellation
    tx_waveform = mapper(tx,mode);
    clear tx;
    
    % noise
    % -----
    noise = 1/sqrt(2)*[randn(1,length(tx_waveform))+sqrt(-1)*randn(1,length(tx_waveform))];
    rx = tx_waveform + 10^(-Es_N0_dB(ii)/20)*noise; % additive white gauss noise
    clear noise;
    clear tx_waveform;
    
    % soft demapper
%     plot(rx,'.');
    llr = demapper(rx,mode);
    clear rx;
       
    % ldpc decode
    H=genH(coderate);
    for nframe = 1:N
        symbol = llr((nframe-1)*9216+1:nframe*9216);
        [recode, iter] = ldpc_decode(symbol,0,H,coderate,col_order);
        for i = 1:9216*coderate
%            outbit(i)=recode(col_order(9216*(1-coderate)+i)+1);
            outbit((nframe-1)*9216*coderate+i)=recode(col_order(9216*(1-coderate)+i)+1);
        end
    clear symbol;
    clear recode;
    end
    nBitErr(ii) = length(find(inbit ~= outbit));
    clear H;
    clear llr;
end

simBer = nBitErr/(N*9216);

semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
axis([0 2 10^-5 1])
grid on
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for QPSK 1/2 coderate')