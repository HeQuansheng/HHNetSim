function Visualization_movieFromRaster(raster,Ntypes,tspan,filename)

%% extract instantanous firing rate
tt=linspace(tspan(1),tspan(end),400);

for i=1:length(raster)
    Unique=sort(unique(raster{i}(:,2)));
    
    for j=1:length(Unique)
        t_spk=raster{i}(raster{i}(:,2)==Unique(j),1)/1000;
        ifr{i}(:,j)= instantfr(t_spk,tt);
    end
    
end

%% Neuron position
temp=randperm(sum(Ntypes));

NN=0;
for i=1:length(Ntypes)
    iset = (1:Ntypes(i))+NN;
    typeidx{i}=temp(iset);
    NN = NN+Ntypes(i);
end
X=[[1:sum(Ntypes)]',rand(sum(Ntypes),1)];
X(typeidx{1},2) = 1*X(typeidx{1},2) + 2;
X(typeidx{3},2) = 3*X(typeidx{3},2) + 0;
X(typeidx{4},2) = 3*X(typeidx{4},2) + 0;

size(X)
length(Ntypes)
length(raster)
Siz=[20,40,20,20];
%%
close all
figure('position',[100 100 800 600],'color','w');

writerObj = VideoWriter(filename);
writerObj.FrameRate =60;
writerObj.Quality=100;
open(writerObj);

for i=1:length(tt)
    
    for j=1:length(raster)
        scatter(X(typeidx{j},1),X(typeidx{j},2), Siz(j),ifr{j}(i,:),'filled');hold on
    end
    h=colorbar;ylabel(h,'Firing rate (Hz)','fontsize',15);set(h,'fontsize',15)
    
    set(gca,'CLim',[0,50]);
    set(gca,'Color',[0.1,0.25,0.4]);
    set(gca,'xtick',[],'ytick',[])
    
    writeVideo(writerObj,getframe(gcf));
    
    drawnow
    
end

close(writerObj);

end

