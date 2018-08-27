function [L,U,P] = ldpc_lu2(A)
% clear all;
% clc;
% A=[  1     1     1     0     0
%      1     0     1     1     0
%      0     0     1     0     1
%      1     0     0     1     0
%      0     0     1     1     1];

[n,n]=size(A);
L=zeros(n,n);
U=zeros(n,n);
P=eye(n);
%piv=zeros(1,n);
AP = A;

for i=1:n
   % AP
    t=find(AP(i:n,i)); % find 1 in column i
    max_col = length(t);
  	min_weight = n + 1; % give a initial value for weight
  	for k = 1:max_col,
  		weight = sum(AP(t(k)+i-1,i:n));
  		if weight < min_weight
  			index = t(k)+i-1;
  			min_weight = weight;
  		end
  	end
  		
  	if( index ~= i)	
  	AP( [i index],:) = AP( [index i],: );
  	P([i index],:) = P([index i],: );
  	L( [i index],:) = L ( [index i],:);
    end
%     AP
  	U(i,i:n) = AP(i,i:n);
    L(i:n,i) = AP(i:n,i);
  	
  	% GAUSS elimination
  	c = find(AP(i+1:n,i));
  	for k = 1 : length(c)
  		c1 = c(k) +i;
  		AP(c1,i:n) = abs ( AP(c1, i:n) - AP(i, i:n));
    end
 
%     for i=1:n
%     AP
%     t=find(AP(i,i:n)); % find 1 in row i
%     max_col = length(t);
%   	min_weight = max_col + 1; % give a initial value for weight
%   	for k = 1:max_col,
%   		weight = sum(AP(i:n, t(k)+i-1));
%   		if weight < min_weight
%   			index = t(k)+i-1;
%   			min_weight = weight;
%   		end
%   	end
%   		
%   	if( index ~= i)	
%   	AP( :,[i index]) = AP( :,[index i] );
%   	P(:,[i index]) = P(:,[index i] );
%     end
%     
%   	U(i,i:n) = AP(i,i:n);
%     L(i:n,i) = AP(i:n,i);
%   	
%   	% GAUSS elimination
%   	c = find(AP(i,i+1:n));
%   	for k = 1 : length(c)
%   		c1 = c(k) +i;
%   		AP(i:n,c1) = abs ( AP(i:n,c1) - AP(i:n,i));
%     end
    
 end    
