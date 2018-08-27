function [bitseq]=demapper(waveform,mode);
%function [bitseq] = demapper(waveform, mode);
%mode: 'bpsk','qpsk','16qam'
%waveform: output data stream
%          amplitude normalize to 1
%bitseq input bit sequence (0, 1)
k_bpsk = 1/sqrt(2);
k_16qam = 1/sqrt(10);

if (mode == 'bpsk')
    for i=1:length(waveform)
        bitseq(i) = real(waveform(i))+imag(waveform(i));
    end
end

if (mode == 'qpsk')
    for i=1:length(waveform)
        bitseq(2*i-1) = real(waveform(i));
        bitseq(2*i) = imag(waveform(i));
    end
end

if (mode == '16qm')
    rx=waveform/k_16qam;
    for i=1:length(rx)
        if(abs(real(rx(i)))<=2)
            bitseq(4*i-3) = real(rx(i));
        elseif(real(rx(i)>2))
            bitseq(4*i-3) = 2*(real(rx(i))-1);
        else
            bitseq(4*i-3) = 2*(real(rx(i))+1);
        end
        bitseq(4*i-1) = abs(real(rx(i))) - 2;
        if(abs(imag(rx(i)))<=2)
            bitseq(4*i-2) = imag(rx(i));
        elseif(imag(rx(i)>2))
            bitseq(4*i-2) = 2*(imag(rx(i))-1);
        else
            bitseq(4*i-2) = 2*(imag(rx(i))+1);
        end
        bitseq(4*i) = abs(imag(rx(i))) - 2;       
    end
end

