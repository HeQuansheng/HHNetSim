function [raster] = GeneratorRaster(tspan,tstim,amp,type,scaleF,N,plot_flag)
if nargin<1;
    N=100;
    dt=0.06;
    tspan=0:dt:1000;
    tstim=[200,700];
    scaleF=0.001;
    amp=[10,0];
    plot_flag=1;
    type='sin';
end

switch type
    case 'external'
        load('modulation.mat');
        Ref=Ch3.values;
        dttemp=Ch3.interval*1000;
        
        tmax=length(Ref)*dttemp;
        tspan=0:dt:tmax;
        
        Ref=downsample(Ref,round(dt./dttemp));
        tspan(end)=[];
        
        spkidx=peakfinder(Ref, 20,-50, 1, 0);
        Ref=spiketriggeredinterpolation(Ref,spkidx,[-0.002,0.005],0.001,1/dt*1000,'linear');
        
        Ref=(Ref-min(Ref))/range(Ref);
        Ref(Ref<0.2)=0;
    case 'sin'
        Ref=GeneratorStimu(tspan,tstim,amp,'sin');
    case 'square'
        Ref=GeneratorStimu(tspan,tstim,amp,'square');

end

Ref=Ref*scaleF;


raster=[];
for i=1:N
    idx= find(Ref>rand(size(Ref)));
    raster=[raster;tspan(idx)',i.*ones(length(idx),1)];
end
raster=sortrows(raster,2);

%%
if plot_flag

    Unique=unique(raster(:,2));
    tt=tspan(1:100:end)./1000;
    for i=1:length(Unique)
        t_spk=raster(raster(:,2)==Unique(i),1);
        ifr(:,i)= instantfr(t_spk./1000,tt);
    end
    
    ifr=mean(ifr,2);

    figure(1),clf
    h1=subplot(211);
    plot(tspan,Ref)
    
    h2=subplot(212);
    [~,H1,~]=plotyy(raster(:,1),raster(:,2),tt.*1000,ifr);
    set(H1,'linestyle','none','marker','.');
    linkaxes([h1,h2],'x')

end


end

