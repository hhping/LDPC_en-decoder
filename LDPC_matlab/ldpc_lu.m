function [L,U,P] = ldpc_lu(A)
%funciton [L,U,P] = ldpc_lu(A)
%A is a NxN matrix, 
%L is lower-triangular matrix, U is a upper-triangular matrix
%P is row-permutation primary matrix
%Do decomposition on A, so that L*U=P*A
%
%                  Editor Chenzy on Dec-15-2009
%Referece artical:
%苏佳宁：适用于数字多媒体广播系统的前向纠错码解码器原理与VLSI实现研究

% clear all;
% clc;
% A=[  1     1     1     0     0
%      1     0     1     1     0
%      0     0     1     0     1
%      1     0     0     1     0
%      0     0     1     1     1];

%Initilizaton
[m,n]=size(A);
if(m~=n)
    error( 'A must be a NxN matrix');
end
U=A;
L=eye(n);
P=eye(n);

%lu decomposition
for i=1:n
%     L
%     U
%     P
%find main cell, minimum row weight
    t = find(U(i:n,i));
%     max_col = length(t);
    min_row_weight = n+1;
    for j=1:length(t)
        row_weight = sum( U(t(j)+i-1,i:n),2);
        if row_weight < min_row_weight
            index = t(j)+i-1;
            min_row_weight = row_weight;
        end
    end
%switch row to get minimal row weight    
    if(index ~= i)
        U([i index],:) = U([index i],:);
        P([i index],:) = P([index i],:);
        if(i>1)
            L([i index],1:i-1)=L([index i],1:i-1);
        end
    end
%Gauss elimination
    c = find(U(i+1:n,i));
    if(length(c)>0)
        for k = 1:length(c)
            c1 = c(k) + i;
            U(c1,i:n) = abs( U(c1,i:n) - U(i,i:n));
            L(c1,i) = abs(L(c1,i) - 1);
        end
    end
end