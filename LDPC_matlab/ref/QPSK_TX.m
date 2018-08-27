%QPSK Transmitter(QPSK_TX.m)
%JC 6/5/05
%Run from editor debug(F5)
%m-file for simulating a QPSK transmitter by modulating with a pseudo
%random bit stream. A serial to parallel conversion of the pseudo random
%bit stream is performed with mapping of two bits per symbol. A cosine and
%sine carrier is configured and the I and Q symbols modulate these
%carriers. The I and Q carriers are combined and time and frequency domain
%plots are provided showing key waveforms at various positions in the QPSK
%transmitter. The simulation uses a serial "passband" approach. Notes and 
%a reference are provided at the end of the m-file.
%===================================================================
clear;
fcarr=2e3;          % Carrier frequency(Hz)
N=10;		        % Number of data bits(bit rate)
fs=16*1e3;		    % Sampling frequency
Fn=fs/2;            % Nyquist frequency
Ts=1/fs;	        % Sampling time = 1/fs
T=1/N;		        % Bit time
randn('state',0);   % Keeps PRBS from changing on reruns
td=[0:Ts:(N*T)-Ts]';% Time vector(data)(transpose)
%===================================================================
% The Transmitter.
%===================================================================
data=sign(randn(N,1))';%transpose
data1=ones(T/Ts,1)*data;
data2=data1(:);

%Serial To Parallel (alternating)
tiq = [0:Ts*2:(N*T)-Ts]';% Time vector for I and Q symbols(transpose)

bs1=data(1:2:length(data));%odd
symbols=ones(T/Ts,1)*bs1;
Isymbols=symbols(:);%I_waveform

bs2=data(2:2:length(data));%even
symbols1=ones(T/Ts,1)*bs2;
Qsymbols=symbols1(:);%Q_waveform

%generate carrier waves
%cosine and sine wave
%2 pi fc t is written as below
twopi_fc_t=(1:fs/2)*2*pi*fcarr/fs; 
a=1;
%phi=45*pi/180
phi=0;
cs_t = a * cos(twopi_fc_t + phi);
sn_t = a * sin(twopi_fc_t + phi);

cs_t=cs_t';%transpose
sn_t=sn_t';
si=cs_t.*Isymbols;
sq=sn_t.*Qsymbols;
sumiq=si+sq;
sumiq=.7*sumiq;%reduce gain to keep output at +/- one
%======================================================================
%Plots
%======================================================================
figure(1)
subplot(3,2,1)
plot(td,data2)
axis([0 1 -2 2]);
grid
xlabel('                                                          Time')
ylabel('Amplitude')
title('Input Data')

subplot(3,2,3)
plot(tiq,Isymbols)
axis([0 1 -2 2]);
grid
xlabel('                                                          Time')
ylabel('Amplitude')
title('I Channel(two bits/symbol) Data')

subplot(3,2,5)
plot(tiq,Qsymbols)
axis([0 1 -2 2]);
grid
xlabel('                                                          Time')
ylabel('Amplitude')
title('Q Channel(two bits/symbol) Data')

subplot(3,2,2)
plot(tiq,si)
axis([.595 .605 -2 2]);
grid
xlabel('                                                          Time')
ylabel('Amplitude')
title('I Channel Modulated Waveform')

subplot(3,2,4)
plot(tiq,sq)
axis([.595 .605 -2 2]);
grid
xlabel('                                                          Time')
ylabel('Amplitude')
title('Q Channel Modulated Waveform')

subplot(3,2,6)
plot(tiq,sumiq)
axis([.595 .605 -2 2]);
grid
xlabel('                                                          Time')
ylabel('Amplitude')
title('QPSK Output Waveform')
%========================================================================
%Take FFT of modulated carrier
%========================================================================
y=sumiq;
NFFY=2.^(ceil(log(length(y))/log(2)));
FFTY=fft(y,NFFY);%pad with zeros
NumUniquePts=ceil((NFFY+1)/2); 
FFTY=FFTY(1:NumUniquePts);
MY=abs(FFTY);
MY=MY*2;
MY(1)=MY(1)/2;
MY(length(MY))=MY(length(MY))/2;
MY=MY/length(y);
f1=(0:NumUniquePts-1)*2*Fn/NFFY;
%=========================================================================
%Plot frequency domain
%=========================================================================
figure(2)
subplot(3,1,1); plot(f1,MY);xlabel('');ylabel('AMPLITUDE');
axis([1500 2500 -.5 2]);%zoom in/out
title('Frequency domain plots');
grid

subplot(3,1,2); plot(f1,20*log10(abs(MY).^2));xlabel('FREQUENCY(Hz)');ylabel('DB');
axis([1500 2500 -100 10]);
grid
title('Modulated QPSK carrier')

%NOTE
%Serial to parallel conversion of a serial bit stream and mapping of
%two bits to a symbol can sometimes be confusing. I will try and explain
%with an example.
%Suppose you have a serial bit stream of ten  0 0 1 1 0 1 1 0 1 1 even # of bits      
                                %odd bits     0   1   0   1   1
                                %even bits      0   1   1   0   1
                                
%The odd bits are the I Channel Data at one half the original serial bit stream
%bit rate. Notice that the amplitudes are +/- one as shown in figure 1.
%The possible combinations are -1 -1, 1 1, -1 1, 1 -1. The amplitudes, in
%theory, should be held at +/- 0.707 to keep the summed output of the 
%QPSK transmitter at a  constant amplitude of +/- one.
%The even bits are the Q Channel Data at one half the original serial bit
%stream bit rate. Same info as above.

%A good reference discussing a QPSK Transmitter and look up tables and Gray
%coding can be found at http://cnx.rice.edu/content/m10042/latest/
