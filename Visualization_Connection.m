function  Visualization_Connection(Ntypes,ConMtx ,pretype,postype,typesHighLight)


% temp=randperm(sum(Ntypes));
% NN=0;
% for i=1:length(Ntypes)
%     iset = (1:Ntypes(i))+NN;
%     typeidx{i}=temp(iset);
%     NN = NN+Ntypes(i);
% end

temp=zeros(1,sum(Ntypes));
temp(3:5:end)=1;
typeidx{1}=find(temp==0);
typeidx{2}=find(temp==1);

MATRIX=zeros(sum(Ntypes),sum(Ntypes));
for i=1:length(Ntypes)
    for j=1:length(Ntypes)
        MATRIX( typeidx{j}, typeidx{i})=ConMtx{i,j};
    end
end

%----------------- define neuronal layout to be viusualize (Start)
MATRIX= sparse(max(MATRIX,MATRIX'));%% be caution that MATRIX is symmetr

%  X = fruchterman_reingold_force_directed_layout(MATRIX);
% [X,~] = kamada_kawai_spring_layout(MATRIX);

X=[[1:size(MATRIX,1)]',rand(size(MATRIX,1),1)];
X(typeidx{2},2) = 1*X(typeidx{2},2) + 2;
% X(typeidx{3},2) = 3*X(typeidx{3},2) + 0;
% X(typeidx{4},2) = 3*X(typeidx{4},2) + 0;

%----------------- define connection to be viusualize (Start)

MATRIX_plot=zeros(size(MATRIX));
MATRIX_plot(typeidx{postype},typeidx{pretype})=MATRIX(typeidx{postype},typeidx{pretype});

%----------------- define connection to be viusualize (End)
figure(1),clf
gplot(MATRIX_plot,X),hold on
set(findobj('Type','line'),'Color',[0.8 0.8 0.8]);

colorma=[1 0 0
                    0 0 1
                    0 1 0
                    0 0 1];
for i=1:length(Ntypes)
    Siz=8;
    if exist('typesHighLight','var')
        if i==typesHighLight
            Siz=15;
        end
    end
    
    plot(X(typeidx{i},1),X(typeidx{i},2),'o','markerfacecolor',colorma(i,:),'markersize',Siz),hold on
end

axis tight
drawnow

%%
figure(2),clf
% MATRIX_shown=[MATRIX(typeidx{1},typeidx{1}),MATRIX(typeidx{1},typeidx{2}),MATRIX(typeidx{1},typeidx{3}),MATRIX(typeidx{1},typeidx{4})
%                               MATRIX(typeidx{2},typeidx{1}),MATRIX(typeidx{2},typeidx{2}),MATRIX(typeidx{2},typeidx{3}),MATRIX(typeidx{2},typeidx{4})
%                               MATRIX(typeidx{3},typeidx{1}),MATRIX(typeidx{3},typeidx{2}),MATRIX(typeidx{3},typeidx{3}),MATRIX(typeidx{3},typeidx{4})
%                               MATRIX(typeidx{4},typeidx{1}),MATRIX(typeidx{4},typeidx{2}),MATRIX(typeidx{4},typeidx{3}),MATRIX(typeidx{4},typeidx{4})];

MATRIX_shown=[MATRIX(typeidx{1},typeidx{1}),MATRIX(typeidx{1},typeidx{2})
                              MATRIX(typeidx{2},typeidx{1}),MATRIX(typeidx{2},typeidx{2})];

imagesc(MATRIX_shown)
colormap(1-gray),hold on

figure(3),clf
NumofTargetReceived=[];
NumofTargetGiveout=[];
for i=1:length(typeidx)
    for j=1:length(typeidx)
        NumofTargetReceived(j,i)=mean(sum(MATRIX(typeidx{j},typeidx{i}),2));
        NumofTargetGiveout(j,i)=mean(sum(MATRIX(typeidx{j},typeidx{i}),1));
    end
end

subplot(121)
imagesc(NumofTargetReceived),hold on
for i=1:size(NumofTargetReceived,1)
    for j=1:size(NumofTargetReceived,2)
        text(i,j,num2str(round(NumofTargetReceived(j,i))),'fontsize',20,'color','r')
    end
end
title('Received Synapse')
xlabel('To')
ylabel('From')
colormap(gray)

subplot(122)
imagesc(NumofTargetGiveout),hold on
for i=1:size(NumofTargetGiveout,1)
    for j=1:size(NumofTargetGiveout,2)
        text(i,j,num2str(round(NumofTargetGiveout(j,i))),'fontsize',20,'color','r')
    end
end
title('Output Synapse')
xlabel('To')
ylabel('From')
colormap(gray)

drawnow
% autoArrangeFigures
end





