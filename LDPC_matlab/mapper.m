function [waveform] = mapper(bitseq, mode)
%function [waveform] = mapper(bitseq, mode);
%bitseq input bit sequence [ 0, 1]
%mode: 'bpsk','qpsk','16qam'
%waveform: output data stream
%          amplitude normalize to 1
k_bpsk = 1/sqrt(2);
k_16qam = 1/sqrt(10);

% y = -2 * x + 1; 
% x = 0 => y = 1;
% x = 1 => y = -1;
if (mode == 'bpsk')
    for i=1:length(bitseq)
        mod = (-2*bitseq(i)+1)*(1+sqrt(-1));
        waveform(i)=k_bpsk*mod;
    end
end

% y = -2 * x + 1; 
if (mode == 'qpsk')
    for i=1:length(bitseq)/2
        Re_symbol = -2*bitseq(2*i-1)+1;
        Im_symbol = -2*bitseq(2*i)+1;
	    mod = Re_symbol + sqrt(-1)*Im_symbol;
	    waveform(i)=k_bpsk*mod;
    end
end

% z = ( - 2 * x + 1) * ( -2 * y + 3);
% x = 0, y = 0 => z = 3;
% x = 0, y = 1 => z = 1;
% x = 1, y = 0 => z = -3
% x = 1, y = 1 => z = -1;
if (mode == '16qm')
    for i=1:length(bitseq)/4
        Re_symbol = (-2*bitseq(4*i-3)+1)*(-2*bitseq(4*i-1)+3);
        Im_symbol = (-2*bitseq(4*i-2)+1)*(-2*bitseq(4*i)+3);
	    mod = Re_symbol + sqrt(-1)*Im_symbol;
        waveform(i)=k_16qam*mod;
    end
end
