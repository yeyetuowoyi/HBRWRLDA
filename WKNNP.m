function [LD_mat_new] = WKNNP( LD_mat, LL_mat, DD_mat, K, a )

[rows,cols]=size(LD_mat);
LD_mat_new=zeros(rows,cols);
y_l=zeros(rows,cols);  
y_d=zeros(rows,cols);  

[sort_L,idx_L] = KNN( LL_mat );  %for miRNA
for i = 1 : rows   
         w=zeros(1,K);
        for j = 1 : K
            w(1,j)=a^(j-1)*sort_L(i,j);
%             w(1,j)=sort_L(i,j)^(j-1);
            y_l(i,:) =  y_l(i,:)+ w(1,j)* LD_mat(idx_L(i,j),:); 
        end                      
            y_l(i,:)=y_l(i,:)/sum(sort_L(i,1:K));  
%             y_l(i,:)=y_l(i,:)/sum(w); 
end

[sort_D,idx_D] = KNN( DD_mat );  %for disease
for i = 1 : cols   
        w=zeros(1,K);

        for j = 1 : K
            w(1,j)=a^(j-1)*sort_D(i,j);
%             w(1,j)=sort_D(i,j)^(j-1);
            y_d(:,i) =  y_d(:,i)+ w(1,j)* LD_mat(:,idx_D(i,j)); 
        end                      
            y_d(:,i)=y_d(:,i)/sum(sort_D(i,1:K)); 
%                 y_d(:,i)=y_d(:,i)/sum(w); 
end

a1=1;
a2=1;
y_ld=(y_l*a1+y_d*a2)/(a1+a2);  

 for i = 1 : rows
     for j = 1 : cols
         LD_mat_new(i,j)=max(LD_mat(i,j),y_ld(i,j));
     end    
 end

end

function [ sort_network ,idx] = KNN( network )
    network= network-diag(diag(network)+1); 
    [sort_network,idx]=sort(network,2,'descend');

end


