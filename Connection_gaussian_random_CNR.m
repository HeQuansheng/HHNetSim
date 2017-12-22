
function MATRIX=Connection_gaussian_random_CNR(Nin,Nout,EdgeN,Nautapse,mode,Lmodel,sigma,CNR_flag,plot_flag);

if nargin<1
    Nin = 1000;
    Nout =1000;
    Lmodel=5000; %length of model um
    sigma=250;
    EdgeN=[20,5];
    Nautapse=10;
    
    CNR_flag=0;
    plot_flag=1;
    mode='gaussian';
end

EdgeN=abs(round(normrnd(EdgeN(1),EdgeN(2),[Nin,1])));

MATRIX=zeros(Nout,Nin);

pos=linspace(0,Lmodel,Nout);
pre=linspace(0,Lmodel,Nin);
for i=1:numel(pre)
    
    dx=pos-pre(i);
    %% using cdf
    %     switch mode
    %         case 'gaussian'
    %             p=normcdf(dx,0,sigma);
    %         case 'uniform'
    %             p=unifcdf(dx,-Lmodel,Lmodel);
    %     end
    %
    %
    %     pp=rand(EdgeN(i));
    %     for j=1:length(pp)
    %         idx=find(p>pp(j),1,'first');
    %         MATRIX(idx,i)=1;% may be 2 contact from one neuron
    %     end
    
    %% using pdf,because of its symmetry
    switch mode
        case 'gaussian'
            p=normpdf(dx,0,sigma);
            p=(p-min(p))./range(p);
        case 'uniform'
            p=ones(size(pos));
    end
    
    Connidx=find(rand(size(pos))<p);
    try
        Connidx=Connidx(randperm(numel(Connidx),EdgeN(i)));
    catch
        warning(['inappropriate fixed Conn. Num, max = ',num2str(numel(Connidx))])
    end
    
    MATRIX(Connidx,i)=1;
end

%%
if  CNR_flag
    pCon=mean(sum(MATRIX,2))./Nout;
    pReci=0.47;% from  Durstewitz 2016;
    k = find(MATRIX>0)';
    % Set probabilities according to the number of neighbours
    % (derived from Perin et al. 2011)
    slope = 20*3.9991/Nout;
    if min(Nout,floor(1/pCon/slope))==0;
        disp('Warning: max_neigh=0 => CNR not applicable!');
    end
    [pair_id, N_neigh] = find_neigh(MATRIX);
    
    % Normalize pCon such that p values to add up to original pCon
    p0 = p_calc(MATRIX, pair_id, N_neigh, pCon, slope);
    Q(:,1:2) = [p0(:,1), p0(:,5)];% CNR 个数，CNR 理论概率
    
    % Get bidirectional pairs
    [pairs, N_neigh_pairs] = find_neigh_recur(MATRIX);
    
    % Select connections randomly according to p (except self-connections)
    pair_id_selected = [];
    for i=1:length(Q(:,1))
        pair_id_act = pair_id(N_neigh == Q(i,1),:);
        pairs_act = pairs(N_neigh_pairs == Q(i,1),:);
        pair_id_old = k(ismember(k, pair_id_act));
        N_rec = floor(Q(i,2)*(pReci/2)*length(pair_id_act));                         % Number of recurrently connected pairs
        N_uni = ceil(Q(i,2)*length(pair_id_act)) - 2*N_rec;                         % Number of unidirectionally connected pairs
        k1 = randperm(length(pairs_act(:,1)));                                      % Randomize pair indices
        
        % Determine number of unidirectional and recurrent connections in old pairs
        dummy = zeros(length(pair_id_old),2);
        dummy(:,1) = ceil(pair_id_old/Nout)';
        dummy(:,2) = pair_id_old'-(dummy(:,1)-1)*Nout;
        pair_id_old_2 = (dummy(:,2)'-1)*Nout + dummy(:,1)'; % counterpart of each pair_id in pair_id_old
        pair_id_old_uni = pair_id_old(~ismember(pair_id_old, pair_id_old_2))';
        pairs_old_rec = [pair_id_old(ismember(pair_id_old, pair_id_old_2))' pair_id_old_2(ismember(pair_id_old, pair_id_old_2))'];
        
        
        % Set bidirectional connections
        if ~isempty(pairs_old_rec) && length(pairs_old_rec(:,1)) > N_rec            % if pairs_old_rec is too long, randomly select from it
            k1_old = randperm(length(pairs_old_rec(:,1)));
            pair_id_selected = [pair_id_selected; pairs_old_rec(k1_old(1:N_rec),1)];
            pair_id_selected = [pair_id_selected; pairs_old_rec(k1_old(1:N_rec),2)];
            
        else                                                                        % if pairs_old_rec is too short, use random pairs from pairs_act
            pair_id_selected = [pair_id_selected; pairs_act(k1(1:N_rec),1)];        % Set bidirectional connections
            pair_id_selected = [pair_id_selected; pairs_act(k1(1:N_rec),2)];        % Set bidirectional connections (counterparts)
        end
        pair_id_act = pair_id_act(~ismember(pair_id_act,pair_id_selected));         % Do not reuse pairs that are already selected
        k1 = randperm(length(pair_id_act));                                         % Randomized pair indices
        
        % Set unidirectional connections
        if length(pair_id_old_uni) > N_uni        % if pair_id_old_uni is too long, randomly select from it
            k1_old = randperm(length(pair_id_old_uni(:,1)));
            pair_id_selected = [pair_id_selected; pair_id_old_uni(k1_old(1:N_uni))];
            
        else                                      % if pair_old is too short, use random pairs from pairs_act
            pair_id_selected = [pair_id_selected; pair_id_act(k1(1:N_uni),:)];
        end;
        
    end
    
    MATRIX(pair_id_selected)=1;
end

%%

if size(MATRIX,1)==size(MATRIX,2)
    MATRIX=MATRIX-diag(diag(MATRIX));
    pair_id_diag = ((1:Nout)-1)*Nout + (1:Nout);
    k1 = randperm(length(pair_id_diag));
    N_diag = round((Nautapse/Nout)*length(pair_id_diag));
    MATRIX(pair_id_diag(k1(1:N_diag)))=1;
end
%%

if plot_flag
    
    figure(1),clf
    imagesc(MATRIX)
    colormap(1-gray)
    caxis([0,1])
    
    sigma=2;
    edges=linspace(-3*sigma,3*sigma,min(size(MATRIX)));
    MATRIX_smooth=zeros(size(MATRIX));
    for i=1:size(MATRIX,2)
        kernel=normpdf(edges,0,sigma).*(edges(2)-edges(1));
        data = conv(MATRIX(:,i),kernel);
        center = ceil(length(edges)/2);
        data=data(center:size(MATRIX,1)+center-1);
        data=(data-min(data))./range(data);
        MATRIX_smooth(:,i)=data(:);
    end
    
    figure(2),clf,hold on
    predist=1000; % pretime
    posdist=1000; % postime
    distIn=linspace(0,Lmodel,Nin);
    distOut=linspace(0,Lmodel,Nout);
    reso=Lmodel./Nout;
    
    for i=1:Nin
        idx=round(max([1,distIn(i)/reso-predist/reso]):min([distIn(i)/reso+posdist/reso,Nout]));
        profile=MATRIX_smooth(idx,i);
        distspan=distOut(idx)-distIn(i);
        plot(distspan,profile,'r')
    end
    
    figure(3),clf
    subplot(211),bar(sum(MATRIX,1)),xlabel('pre'),title(['output #  ',num2str(mean(sum(MATRIX,1)))])
    subplot(212),bar(sum(MATRIX,2)),xlabel('pos'),title(['Receive # ',num2str(mean(sum(MATRIX,2)))])
    
    if CNR_flag
        % Analyse connectivity
        pCon=mean(sum(MATRIX,2))./Nout;
        [pair_id_out, N_neigh_out] = find_neigh(MATRIX);
        p = p_calc(MATRIX, pair_id_out, N_neigh_out, pCon, slope);
        
        % Compute fraction of reciprocal connections
        test = (MATRIX>0) + (MATRIX>0)';
        pReci = length(test(test>1)) / length(MATRIX(MATRIX>0))
        
        figure(4),clf
        plot(Q(:,1),Q(:,2),'b-o'),hold on
        plot(p(:,1),p(:,5),'r-s')
        xlabel('# of common neighber')
        ylabel('P_{conn}')
        legend('before','after CNR');
        
    end
    
end

% save('PARAMS_connection.mat','ConMtx','PCidx','FSidx','LTSidx','MATRIX')

end


% -------------------------------------------------------------------------
% ---------------------  Auxillary functions  -----------------------------
% -------------------------------------------------------------------------

function [pair_id, N_neigh] = find_neigh(X)
% function to compute the number of common neighbours for each pair of
% neurons connected by X
% 输出：都是向量
%

Nin = length(X(1,:));
Nout = length(X(:,1));

% Take out self-connections
ind_diag = ((1:Nout)-1)*Nout + (1:Nout);
X(ind_diag) = zeros(1,Nout);

% Find pairs which share at least one common neighbour
X2 = (X + X')>0;             % make matrix symmetrical, as only pairs are important, not directionality of connections
dummy = cumsum(X2>0,2);%作为突触后每个接受多少突触前
pairs = zeros(1,2);
kk=0;
for i=2:max(max(dummy))    % loop up to maximal number of neighbours of any of the cells
    idx = find(dummy(:,end) == i);    % neurons with i outgoing connections
    neighbours = zeros(length(idx),i);
    for j=1:length(idx)
        neighbours(j,:) = find(X2(idx(j),:)>0);
        pairs(kk+1:kk+nchoosek(length(neighbours(j,:)),2),:) = nchoosek(neighbours(j,:),2);     % Get all combinations of neurons that project to that neuron as pairs
        % 有意思的函数 nchoosek
        kk = kk+nchoosek(length(neighbours(j,:)),2);
    end;
end;

% Make indices for bidirectional connections for each pair
% and compute number of neighbours
pair_id = [(pairs(:,1)-1)*Nout + pairs(:,2); (pairs(:,2)-1)*Nout + pairs(:,1)];
[pair_id, ~, pair_ind_redund] = unique(sort(pair_id));
N_neigh = accumarray(pair_ind_redund,1);

% Merge with random pairs, keep only unique pairs (as neighboured)
out_ind = ~ismember(1:Nin*Nout, pair_id);
N_neigh = [zeros(length(find(out_ind)),1); N_neigh];
[pair_id, sort_idx] = sort([find(out_ind)'; pair_id]);
N_neigh = N_neigh(sort_idx);

% Set self-connections to "-1 neighbours"
N_neigh(ind_diag) = -1;


end

function [pair_id, N_neigh] = find_neigh_recur(X)
% function to compute the number of common neighbours for each pair of
% neurons connected by X

Nout = length(X(:,1));

% Take out self-connections
ind_diag = ((1:Nout)-1)*Nout + (1:Nout);
X(ind_diag) = zeros(1,Nout);

% Find pairs which share at least one common neighbour
X2 = (X + X')>0;                            % make matrix symmetrical, as only pairs are important, not directionality of connections
dummy = cumsum(X2>0,2);
pairs = zeros(1,2);
kk=0;
for i=2:max(max(dummy))
    idx = find(dummy(:,end) == i);          % neurons with i outgoing connections
    neighbours = zeros(length(idx),i);
    for j=1:length(idx)
        neighbours(j,:) = find(X2(idx(j),:)>0);
        pairs(kk+1:kk+nchoosek(length(neighbours(j,:)),2),:) = nchoosek(neighbours(j,:),2);
        kk = kk+nchoosek(length(neighbours(j,:)),2);
    end;
end;

% Make indices for bidirectional connections for each pair
% and compute number of neighbours
pair_id(:,1) = (pairs(:,1)-1)*Nout + pairs(:,2);
pair_id(:,2) = (pairs(:,2)-1)*Nout + pairs(:,1);
[~, pair_ind_unique, pair_ind_redund] = unique(pair_id(:,1));
pair_id = pair_id(pair_ind_unique,:);
N_neigh = accumarray(pair_ind_redund,1);

% Get connections on the lower triangluar part of X
X_tril = tril(ones(size(X)));
tril_ind = find((X_tril>0));
out_ind = ~ismember(tril_ind, pair_id(:,1));
pair_id_zero(:,1) = tril_ind(out_ind);

% Reconstruct connections on the upper triangluar part of X
pairs_zero(:,1) = ceil(pair_id_zero(:,1)/Nout);
pairs_zero(:,2) = pair_id_zero(:,1)-(pairs_zero(:,1)-1)*Nout;
pair_id_zero(:,2) = (pairs_zero(:,2)-1)*Nout + pairs_zero(:,1);

% Add zero neigbour pairs (all the rest of X_tril)
N_neigh = [zeros(length(find(out_ind)),1); N_neigh];
pair_id = [pair_id_zero; pair_id];
[~, sort_idx] = sort(pair_id(:,1));
pair_id = pair_id(sort_idx,:);
N_neigh = N_neigh(sort_idx);

% Set self-connections to "-1 neighbours"
N_neigh(ismember(pair_id(:,1),ind_diag)) = -1;

end

function p = p_calc(X, pair_id, N_neigh, pCon, slope)
% function to analyse connectivity
% p(:,1) ---- Num of common neighbor,NC
% p(:,2) ---- 含有对应NC的细胞对
% p(:,3) ---- 含有对应NC，且有连接的细胞对
% p(:,4) ---- 2/3, 对应NC的实际概率.
% p(:,5) ---- 对应NC的理论概率.

Nout = length(X(:,1));

% Take out self-connections
ind_diag = ((1:Nout)-1)*Nout + (1:Nout);
X(ind_diag) = zeros(1,Nout);

p(:,1) = unique(N_neigh(N_neigh>=0));    % do not use self-connections
pair_id_selected = find(X>0);
N_neigh_selected = N_neigh(ismember(pair_id, pair_id_selected));% 分析有连接的pair
for i=1:length(p(:,1))
    pair_act = pair_id(N_neigh == p(i,1));
    p(i,2) = length(pair_act);
    p(i,3) = sum(N_neigh_selected==p(i,1));
    p(i,4) = sum(N_neigh_selected==p(i,1)) / length(pair_act);
end;

% Normalize pCon such that p values to add up to original pCon
off = fminsearch(@(N0) p_min(N0, p, pCon, slope), max(p(:,1)));
N1 = floor(min(max(p(:,1)),  off+1/slope/pCon));
N0 = max(min(p(:,1)), ceil(off));
ind = ismember(p(:,1), N0:N1);
p(:,5) = zeros(size(p(:,1)));
p(ind,5) = pCon*slope*(p(ind,1)-off);
p(find(ind,1,'last')+1:end,5) = 1;
if ~isempty(find(p(:,5)<0, 1))
    disp('Warning: negative probabilities!');
    n_idx=find(p(:,5)<0);
    disp(['i=' num2str(n_idx)]);
elseif ~isempty(find(p(:,5)>1, 1))
    disp('Warning: probabilities > 1');
    p_idx=find(p(:,5)<0);
    disp(['i=' num2str(p_idx)]);
end

end

function res = p_min(N0, p, pCon, slope)
% function to be minimized to compute offset for p_Neighbours

N1 = round(min(max(p(:,1)),  N0+1/slope/pCon));

N_2 = p(p(:,1)>N0 & p(:,1)<=N1,1);
pN_2 = p(p(:,1)>N0 & p(:,1)<=N1,2)/sum(p(:,2));
pN_3 = p(p(:,1)>N1,2)/sum(p(:,2));

res = abs(sum([pN_2*slope.*(N_2-N0); pN_3]) -1);


end
% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg
% and BCCN Heidelberg-Mannheim