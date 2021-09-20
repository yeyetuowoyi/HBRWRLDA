function [HGP] = Highgraph( Y ,S )
[n,m]=size(Y);
DV=diag(sum(Y,2));
WS=sum(S,2);
W=zeros(n,m);

for i=1:n
    for j=1:m
        W(i,j)=Y(i,j)*WS(i);
%          W(i,j)=H(i,j)*WS(i);
% W(i,j)=WS(i);
    end
end
DVe=diag(sum(W));
We=diag(ones(m,1)/m);

HGP=DV^-1*Y*We*DVe^-1*W';
end