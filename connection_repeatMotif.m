function [ConMtx,typeidx] = connection_repeatMotif(NtypesMotif,NumCopy,Ntypes,pCon,pConMotif,wgt)
if nargin<1
    
%     % testting data
%     NtypesMotif=[5,2,2,0];
%     Ntypes= [0,0,0,2];
%     NumCopy=10;
%     
%     pConMotif=ones(length(NtypesMotif),length(NtypesMotif)).*1;
%     pCon=ones(length(NtypesMotif),length(NtypesMotif)).*0.02;
%     % cross-column/scatter connection
%     pCon(1,2)=1;
%     pCon(4,3)=1;
%     
%     % within column
%     pConMotif(1,1)=20./NtypesMotif(1);
%     pConMotif(3,3)=1./NtypesMotif(3);
%     


    NtypesMotif=[50,80,10,0];
    Ntypes= [0,0,0,20];
    NumCopy=10;
    
    pConMotif=ones(length(NtypesMotif),length(NtypesMotif)).*0.02;
    pCon=ones(length(NtypesMotif),length(NtypesMotif)).*0.02;
    % cross-column/scatter connection
    pCon(1,1)=0.1*Ntypes(1)./Ntypes(1);
    pCon(1,2)=0.1*Ntypes(2)./Ntypes(1);
    pCon(1,3)=0.1*Ntypes(3)./Ntypes(1);
    pCon(1,4)=0.1*Ntypes(4)./Ntypes(1);
    
    pCon(2,1)=0.2*Ntypes(1)./Ntypes(2);
    pCon(2,2)=0.2*Ntypes(2)./Ntypes(2);
    pCon(2,3)=0.1*Ntypes(3)./Ntypes(2);
    pCon(2,4)=0.1*Ntypes(4)./Ntypes(2);
    
    pCon(3,1)=0.1*Ntypes(1)./Ntypes(3);
    pCon(3,2)=0.1*Ntypes(2)./Ntypes(3);
    pCon(3,3)=0.0*Ntypes(3)./Ntypes(3);
    pCon(3,4)=0.8*Ntypes(4)./Ntypes(3);
    
    pCon(4,1)=0.01*Ntypes(1)./Ntypes(4);
    pCon(4,2)=0.01*Ntypes(2)./Ntypes(4);
    pCon(4,3)=0.1*Ntypes(3)./Ntypes(4);
    pCon(4,4)=0.1*Ntypes(4)./Ntypes(4);
    
    % within column
    pConMotif(1,1)=20./NtypesMotif(1);
    pConMotif(1,2)=20./NtypesMotif(1);
    pConMotif(1,3)=10./NtypesMotif(1);
    
    pConMotif(2,1)=10./NtypesMotif(2);
    pConMotif(2,2)=20./NtypesMotif(2);
    pConMotif(2,3)=10./NtypesMotif(2);
    
    pConMotif(3,1)=1./NtypesMotif(3);
    pConMotif(3,2)=1./NtypesMotif(3);
    pConMotif(3,3)=1./NtypesMotif(3);
    
    
    wgt(1,1)={[log([0.001,0.5]),0.1,0.5]};
    wgt(1,2)={[log([0.001,0.5]),0.1,0.5]};
    wgt(1,3)={[log([0.0001,0.5]),2,3]};
    wgt(1,4)={[log([0.0001,0.5]),2,3]};
    
    wgt(2,1)={[log([0.001,0.5]),0.1,0.5]};
    wgt(2,2)={[log([0.001,0.5]),0.1,0.5]};
    wgt(2,3)={[log([0.0001,0.5]),2,3]};
    wgt(2,4)={[log([0.0001,0.5]),2,3]};
    
    wgt(3,1)={[log([50,0.5]),100,200]};
    wgt(3,2)={[log([50,0.5]),100,200]};
    wgt(3,3)={[log([50,0.5]),100,200]};
    wgt(3,4)={[log([50,0.5]),100,200]};
    
    wgt(4,1)={[log([50,0.5]),100,200]};
    wgt(4,2)={[log([50,0.5]),100,200]};
    wgt(4,3)={[log([50,0.5]),100,200]};
    wgt(4,4)={[log([50,0.5]),100,200]};

end

%% make BaseMatrix contain repeat motif
motifidx=find(NtypesMotif~=0);
for i=1:length(motifidx)
    for j=1:length(motifidx)
        subConMtx{j,i}=rand(NtypesMotif(motifidx(j)),NtypesMotif(motifidx(i)))>(1-pConMotif(motifidx(i),motifidx(j)));
    end
end
Motif=cell2mat(subConMtx);
BaseMATRIX = kron(eye(NumCopy),Motif);
%% insert BaseMatrix into MATRIX
NTYPES=NtypesMotif.*NumCopy+Ntypes;
NtypesSimple=NtypesMotif+Ntypes;
MATRIX=zeros(sum(NTYPES),sum(NTYPES));

NN=0;
for i=1:length(Ntypes)
    if NtypesMotif(i)>0
        for j=1:NumCopy
            iset = (1:NtypesMotif(i))+NN;
            increa=NtypesMotif(i)*(j-1);
            iset=iset+increa;
            typeidx{i}{j}=sort(iset(:));
        end
    else
        iset = (1:Ntypes(i))+NN;
        typeidx{i}{1}=sort(iset(:));
        for qq=2:NumCopy
            typeidx{i}{qq}=[];
        end
    end
    
    NN=NN+NTYPES(i);
end

NN=0;
for i=1:length(NtypesMotif)
    for j=1:NumCopy
        iset = (1:NtypesMotif(i))+NN;
        increa=sum(NtypesMotif)*(j-1);
        iset=iset+increa;
        typeidx_Base{i}{j}=sort(iset(:));
    end
    NN = NN+NtypesMotif(i);
end


for i=1:length(NtypesSimple)
    for j=1:length(NtypesSimple)
        for k=1:NumCopy
            for q=1:NumCopy
                
                if ~isempty(typeidx_Base{i}{k})&&~isempty(typeidx_Base{j}{q})&&q==k
                    MATRIX(typeidx{j}{k},typeidx{i}{q})=BaseMATRIX(typeidx_Base{j}{k},typeidx_Base{i}{q});
                else
                    MATRIX(typeidx{j}{k},typeidx{i}{q})=rand(numel(typeidx{j}{k}),numel(typeidx{i}{q}))>(1-pCon(i,j));
                end
                
            end
        end
    end
end

MATRIX=MATRIX-diag(diag(MATRIX));

%%

for i=1:length(NTYPES)
    for j=1:length(NTYPES)
        X=zeros(NTYPES(i),NTYPES(i));
        idxSource=[];
        idxTarget=[];
        for k=1:NumCopy
            idxSource=[idxSource;typeidx{i}{k}];
            idxTarget=[idxTarget;typeidx{j}{k}];
        end
        X= MATRIX(idxTarget,idxSource);
        X(X>0)=GeneratorRand(sum(X(:)), wgt{i,j}(1), wgt{i,j}(2), wgt{i,j}(3), wgt{i,j}(4), 1);
        
        ConMtx(i,j)={X};
        
    end
end

% save('connection.mat','ConMtx','NTYPES','typeidx')

%%
figure(1),clf
imagesc(MATRIX),hold on
colormap(1-gray)
NN=0;
for i=1:length(NTYPES)-1
    plot([NTYPES(i);NTYPES(i)]+NN,[0;sum(NTYPES)],'r--'),hold on
    plot([0;sum(NTYPES)],[NTYPES(i);NTYPES(i)]+NN,'r--')
    NN=NN+NTYPES(i);
end
axis tight

% figure(2),clf
% imagesc(BaseMatrix)
% colormap(1-gray)
end

