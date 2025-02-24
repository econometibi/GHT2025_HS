function Wn = normw(Wn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Wn = Wn - diag(diag(Wn)); % ensure diagonal is zero
% sums = sum(Wn');
D = inv(diag(sum(Wn')));

Wnorm = D*Wn; % row normalize W matrix
Wn=Wnorm;
end

