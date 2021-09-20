lncSim = load ('.\Dataset\lncRNAsimilarity.txt');
interaction = importdata ('.\Dataset\known_lncRNA_disease_interaction.txt');
disSim = load('.\Dataset\diseasesimilarity.txt');
%parameter
k=40;
w=0.8;
L1=2;
L2=2;
alpha=0.5;


A_ori=interaction;

[nl,nd]=size(A_ori);


A= WKNNP( A_ori,  lncSim,disSim, k, w);

[P1] = Highgraph( A ,lncSim );
[P2] = Highgraph( A' ,disSim );


[Rt]=BR(A,P1,P2,L1,L2,alpha);

%loocv
index=find(A_ori==1);
for u=1:length(index)
    
    Y=A_ori ;
    Y(index(u))=0;


    LD= WKNNP( Y,   lncSim,disSim,  k, w);
    [P1] = Highgraph( LD ,lncSim );
    [P2] = Highgraph( LD' ,disSim );

    [Score]=BR(LD,P1,P2,L1,L2,alpha);
    Rt(index(u))=Score(index(u));
end
    pre_label_score = Rt(:);
    label_y = A_ori(:);
    auc= roc_1(pre_label_score,label_y,'red');
