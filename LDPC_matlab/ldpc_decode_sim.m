%function [R iter]=ldpc_decode_fix(llr,SNR,H,coderate,col_order)
%function [R iter]=ldpc_decode(llr,SNR,H,coderate,col_order);
%llr : log-like ratio
%SNR : SNR value estimation
%H   : Check Matrix
%col_order : Reranged Column Order
%R   : decoder result
%ber : bit error rate

% llr width
% llr_width = 5;
load sim.mat
load llr.mat

iter_num = 20;
alpha = 0.625;
ber = 0;
if coderate == 0.5
    row_weight = 6;
    row_num = 18;
elseif coderate == 0.75
    row_weight = 12;
    row_num = 9;
end

% store first 18 row only for 0.5 coderate
% store first 9 row for 0.75 coderate

% initilize ram index and offset value
% 18x6 store index and offset
if coderate == 0.5
    t=zeros(18,6);
    for i=1:18
        t(i,:)=find(H(i,:))';
    end
    index=mod(t,36);
    offset=fix(t/36)+1;
    % change 0 to 36
    n=find(index==0);
    index(n)=36;
    offset(n)=offset(n)-1;
% 9x12 store index and offset    
elseif coderate == 0.75
    t=zeros(9,12);
    for i=1:9;
        t(i,:)=find(H(i,:))';
    end
    index=mod(t,36);
    offset=fix(t/36)+1;
    % change 0 to 36
    n=find(index==0);
    index(n)=36;
    offset(n)=offset(n)-1;
end

% initilize seq_ram
for i =1:256
    for j =1:36
        seq_ram(j,i) = llr(36*(i-1)+j);
%         temp = floor(llr(36*(i-1)+j)*2^llr_width;
%	 seq_ram(j,i) = temp;
    end
end

% H matrix store, t(i,j) non-zero t_group no-zero m group 
% for i = 1:256
%     for j =1:256
%         temp = H((18*i-17):(18*i),(36*j-35):(36*j));
%         t(i,j)=nnz(temp);
%     end
% end
% t_group = find(t~=0);

compresslr=zeros(row_num,256,row_weight+4);

for iter = 1:iter_num
% Horizontal step: collect lq information    
    for i = 1:256
        for j = 1:row_num
            for k = 1:row_weight
                x = index(j,k);
                y = mod(offset(j,k)+i-1,256);
                if y == 0
                    y = 256;
                end
                lqij(k) = seq_ram(x,y);
%                 if x == 36 & y == 6
%                    x
%                    y
%                    seq_ram(x,y)
%                 end
            end    
            
            % Decompress : generate lr information
            for k = 1:row_weight
                if k == compresslr(j,i,3)
                    temp_min = compresslr(j,i,2);
                else
                    temp_min = compresslr(j,i,1);
                end
                sign_temp = mod( compresslr(j,i,k+4)+compresslr(j,i,4), 2);
                if sign_temp == 1
                    lrij(k) = -1 * temp_min;
                else
                    lrij(k) = temp_min;
                end
            end
            
            % Retrive lqij' information
            if iter == 1
                lqij = lqij;
            else
                lqij = lqij -lrij;
            end
            
            % Min-sum function
            % Sign of lqij
            sign_xor = 0;
            for k = 1:row_weight
                if lqij(k) < 0;
                    sign_xor = mod (sign_xor + 1,2);
                    sign_lq(k) = 1;
                else
                    sign_lq(k) = 0;
                end
            end
            
            % Minimum value
            if abs(lqij(1)) < abs(lqij(2))
                min_lq = abs(lqij(1));
                less_lq = abs(lqij(2));
                loc_lq = 1;
            else
                min_lq = abs(lqij(2));
                less_lq = abs(lqij(1));
                loc_lq = 2;
            end
            for k = 3:row_weight
                if abs(lqij(k)) < min_lq
                    less_lq = min_lq;
                    min_lq = abs(lqij(k));
                    loc_lq = k;
                elseif abs(lqij(k)) < less_lq
                    less_lq = abs(lqij(k));
                end
            end
            min_lq = round(alpha * min_lq);
            less_lq = round(alpha * less_lq);
            compresslr(j,i,1) = min_lq;
            compresslr(j,i,2) = less_lq;
            compresslr(j,i,3) = loc_lq;
            compresslr(j,i,4) = sign_xor;
            compresslr(j,i,5:end) = sign_lq(1:row_weight);
                 
             % Decompress : generate lr information
            for k = 1:row_weight
                if k == compresslr(j,i,3)
                    temp_min = compresslr(j,i,2);
                else
                    temp_min = compresslr(j,i,1);
                end
                sign_temp = mod( compresslr(j,i,k+4)+compresslr(j,i,4), 2);
                if sign_temp == 1
                    lrij(k) = -1 * temp_min;
                else
                    lrij(k) = temp_min;
                end
            end
            
            %add new lr back to lq, then write lq back to ram_seq
            lqij = lqij + lrij;
            
            for k = 1:row_weight
                x = index(j,k);
                y = mod(offset(j,k)+i-1,256);
                if y == 0
                    y = 256;
                end
                seq_ram(x,y) = lqij(k);
            end
        end
    end
    for i =1:256
        for j =1:36
            llr(36*(i-1)+j) = seq_ram(j,i) ;
        end
    end
% Hard decision output data
    for i =1:length(llr)
        R(i) = (llr(i) < 0);
    end
    PC = mod((H*R'),2);
    if nnz(PC) == 0
        break;
    end
end

    



