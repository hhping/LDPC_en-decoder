%function [L,U,P] = ldpc_lu(A)
clear all;
clc;
A=[  1     1     1     0     0
     1     0     1     1     0
     0     0     1     0     1
     1     0     0     1     0
     0     0     1     1     1];

%initialization
% [n,n]=size(A);
% L=zeros(n,n);
% U=zeros(n,n);
% P=eye(n);
% %piv=zeros(1,n);
% AP = A;

%lu decompose
% for i=1:n
%    % AP
%     t=find(AP(i:n,i)); % find 1 in column i
%     max_col = length(t);
%   	min_weight = n + 1; % give a initial value for weight
%   	for k = 1:max_col,
%   		weight = sum(AP(t(k)+i-1,i:n));
%   		if weight < min_weight
%   			index = t(k)+i-1;
%   			min_weight = weight;
%   		end
%   	end
%   		
%   	if( index ~= i)	
%   	AP( [i index],:) = AP( [index i],: );
%   	P([i index],:) = P([index i],: );
%   	L( [i index],:) = L ( [index i],:);
%     end
% %     AP
%   	U(i,i:n) = AP(i,i:n);
%     L(i:n,i) = AP(i:n,i);
  	
%   	% GAUSS elimination
%   	c = find(AP(i+1:n,i));
%   	for k = 1 : length(c)
%   		c1 = c(k) +i;
%   		AP(c1,i:n) = abs ( AP(c1, i:n) - AP(i, i:n));
%     end
 
[m,n]=size(A);

U=A;
L=eye(m);
P=eye(m);

for i=1:m
    i
    L
    U
    P
    %find main cell
    t =(find(U(i:end,i)==1));
    max_col = length(t);
    min_row_weight1 = n+1;
    for k =1:max_col
        row_weight = sum( U(t(k)+i-1,i:end),2);   
        if  row_weight < min_row_weight1
            row_index = t(k)+i-1;
            min_row_weight1 = row_weight;
        end
     end
     
     if(row_index ~= i)
     U([i row_index],:) =  U([row_index i],:);
     P([i,row_index],:) =  P([row_index i],:);    
     if(i>1)
          L([i row_index],1:i-1)=L([row_index i],1:i-1);
         end
     end
     
     %Gauss elimination           
     y2=(find(U(i+1:end,i)==1));
     k2 = length(y2);
     if k2>0
         for k3=1:k2
               k4=y2(k3)+i;
               U(k4,i:end)=abs( U(k4,i:end) - U(i,i:end));
               L(k4,i) = abs(L(k4,i) -1);
         end  
      end
 end