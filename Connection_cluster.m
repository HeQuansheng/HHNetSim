function [ConMtx] = Connection_cluster(Ntypes,pCon,pSelfCon,pRec,wgt,cluster_flag,distr_flag)

for i=1:length(Ntypes)
    for j=1:length(Ntypes)
        
        switch cluster_flag(i,j)
            case 1 % CNR
                [MATRIX,~] = CommonNeighbour_Recur(0, Ntypes(j), Ntypes(i), pCon(i,j), pSelfCon(i,j), pRec(i,j),1);
                
            case 0 % randomly
                MATRIX = Connection_Rand(Ntypes(j), Ntypes(i), pCon(i,j));
                
            case 2 % gaussian
                MATRIX = Connection_gaussian(Ntypes(i),Ntypes(j),5000,250,1);
        end;
        MATRIX=MATRIX>0;
        MATRIX(MATRIX>0)=GeneratorRand(sum(MATRIX(:)), wgt{i,j}(1), wgt{i,j}(2), wgt{i,j}(3), wgt{i,j}(4), distr_flag);
       
        ConMtx(i,j)={MATRIX};
        
    end
    
end

% save('connection.mat','ConMtx','Ntypes')

end


function [X]=Connection_Rand(Nin,Nout,pCon)
% k1 = randperm(Nin*Nout);
% k = k1(1:round(pCon*Nin*Nout));
% X=zeros(Nin,Nout);
% idc=Nsyn+(1:length(k));
% X(k)=idc;
X=rand(Nin,Nout)>(1-pCon);


end


function MATRIX = Connection_gaussian(Nin,Nout,Lmodel,sigma,sametype_flag,scalefactor)
% 1--L2PC,2--L5PC,3--LTS, 4--FS

MATRIX=zeros(Nout,Nin);

pos=linspace(0,Lmodel,Nout);
pre=linspace(0,Lmodel,Nin);
for i=1:numel(pre)
    % gaussian distance
    dx=pos-pre(i);
    p=normpdf(dx,0,sigma);
    con=rand(size(pos)).*scalefactor <p;
    MATRIX(:,i)=con;
end

if sametype_flag
    MATRIX=MATRIX-diag(diag(MATRIX));
end

end


