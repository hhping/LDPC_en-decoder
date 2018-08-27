function [H]=genH(coderate);
%function [H]=genH(coderate);
%coderate can be set 1/2 (0.5) or 3/4 (0.75)
%when coderate equals 0.5
%  H would be a matrix 4608x9216
%when coderate equlas 0.75
%  H must be a matrix 2304x9216
%
%                  Editor Chenzy on Dec-15-2009

if(coderate == 0.5)
    rows = 4608;
    cols = 9216;
%row_flag(1:rows)=0;
    HT=zeros(rows,cols);

%initialization    
    aij = zeros(18,9216);
    r = zeros(18,6);
%look up table from protocol    
    r(1,:)= [0 6 12 18 25 30];
    r(2,:)= [0 7 19 26 31 5664];
    r(3,:)= [0 8 13 20 32 8270];
    r(4,:)= [1 6 14 21 3085 8959];
    r(5,:)= [1 15 27 33 9128 9188];
    r(6,:)= [1 9 16 34 8485 9093];
    r(7,:)= [2 6 28 35 4156 7760];
    r(8,:)= [2 10 17 7335 7545 9138];
    r(9,:)= [2 11 22 5278 8728 8962];
    r(10,:)= [3 7 2510 4765 8637 8875];
    r(11,:)=[3 4653 4744 7541 9175 9198];
    r(12,:)=[3 23 2349 9012 9107 9168];
    r(13,:)=[4 7 29 5921 7774 8946];
    r(14,:)=[4 7224 8074 8339 8725 9212]; 
    r(15,:)=[4 4169 8650 8780 9023 9159];
    r(16,:)=[5 8 6638 8986 9064 9210];
    r(17,:)=[5 2107 7787 8655 9141 9171];
    r(18,:)=[5 24 5939 8507 8906 9173];
 %child-matrix gen 18x9216
    j = 1;
    for k=1:18
	    for i=1:cols
  	        if  (i-1) == r(k,j)
    	        aij(k,i) = 1;
    	        if j<6
                    j = j + 1; 
                end
             end
        end
        j = 1;
    end    
%QC-shift generate full matrix
% A1   A2 A3 ... A256
% A256 A1 A2     A255
% A2   A3 A4     A1
    for i = 0:255
        for j = 0:255
            HT( (i*18+1):(i*18+18), (j*36+1):(j*36+36))= aij( :,( mod(j+256-i,256)*36+1):( mod(j+256-i,256)*36+36) );
        end
    end
   H=HT;
end

if(coderate == 0.75)
    rows = 2304;
    cols = 9216;
%row_flag(1:rows)=0;
    HT=zeros(rows,cols);

%initialization    
    aij = zeros(18,9216);
    r = zeros(9,12);
%look up table from protocol    
    r(1,:)= [0 3 6 12 16 18 21 24 27 31 34 7494];
    r(2,:)= [0 4 10 13 25 28 5233 6498 7018 8358 8805 9211];
    r(3,:)= [0 7 11 19 22 6729 6831 7913 8944 9013 9133 9184];
    r(4,:)= [1 3 8 14 17 20 29 32 5000 5985 7189 7906];
    r(5,:)= [1 9 4612 5523 6456 7879 8487 8952 9081 9129 9164 9214];
    r(6,:)= [1 5 23 26 33 35 7135 8525 8983 9015 9048 9154];
    r(7,:)= [2 3 30 3652 4067 5123 7808 7838 8231 8474 8791 9162];
    r(8,:)= [2 35 3774 4310 6827 6917 8264 8416 8542 8834 9044 9089];
    r(9,:)= [2 15 631 1077 6256 7859 8069 8160 8657 8958 9094 9116];
 %child-matrix gen 9x9216
    j = 1;
    for k=1:9
	    for i=1:cols
  	        if  (i-1) == r(k,j)
    	        aij(k,i) = 1;
    	        if j<12
                    j = j + 1; 
                end
             end
        end
        j = 1;
    end    
%QC-shift generate full matrix
% A1   A2 A3 ... A256
% A256 A1 A2     A255
% A2   A3 A4     A1
    for i = 0:255
        for j = 0:255
            HT( (i*9+1):(i*9+9), (j*36+1):(j*36+36))= aij(1:9,( mod(j+256-i,256)*36+1):( mod(j+256-i,256)*36+36) );
        end
    end
   H=HT;
end