function y=GeneratorStimu(tspan,tstim,amp,type,intensity,N)
% function used to generate arbitory waveform at given time period
% 
% INPUTS:
% tstim        : time point to apply stimulus, [on,off]------- n x 2 matrix
% amp            : stimulus ampitude,                     [A1..Arest]----- (n+1) x 1 vector
% intensity: in Hz, then tspan should be in second (better keep in ms in future version)
% 
% OUPUTS:
% y                  : a stimulus vector same size with tspan
% 
% EXAMPLE:
% tspan=linspace(0,1000,10000); 
% tstim=[100,400;400,800;800,tspan(end)];
% amp=[-120,40,-80,-90];
% y=GeneratorStimu(tspan,tstim,amp,'square');
% plot(tspan,y)
%
% He Quansheng
% quanshe@mail.bnu.edu.cn

dt=tspan(2)-tspan(1);

if nargin==4
     intensity=zeros(size(tstim,1),1)+100;
     N=10;
end


switch type
    case 'square'
        y=zeros(size(tspan))+amp(end);% resting state
        for i=1:size(tstim,1)
            stimidx=round(tstim(i,1)/dt:tstim(i,2)/dt);
            stimidx(stimidx>length(tspan))=[];
            y(stimidx)=amp(i);
        end
        
    case 'sin'
        y=zeros(size(tspan))+amp(end);% resting state
        for i=1:size(tstim,1)
            stimidx=round(tstim(i,1)/dt:tstim(i,2)/dt);
            stimidx(stimidx>length(tspan))=[];
            freq=1./(tstim(i,2)-tstim(i,1))./2;
            y(stimidx)=amp(i).*sin(2*pi*freq.*(tspan(stimidx)-tspan(stimidx(1))));
        end
        
   case 'rand'
       y=zeros(size(tspan))+amp(end);% resting state
        for i=1:size(tstim,1)
            stimidx=round(tstim(i,1)/dt:tstim(i,2)/dt);
            stimidx(stimidx>length(tspan))=[];
            y(stimidx)=amp(i).*rand(size(stimidx));
        end
        
    case 'ramp'
        y=zeros(size(tspan))+amp(end);% resting state
        for i=1:size(tstim,1)
            stimidx=round(tstim(i,1)/dt:tstim(i,2)/dt);
            stimidx(stimidx>length(tspan))=[];
            y(stimidx)=amp(i).*linspace(0,1,numel(stimidx));
        end
        
    case 'chirp'
        y=zeros(size(tspan))+amp(end);% resting state
        for i=1:size(tstim,1)
            stimidx=round(tstim(i,1)/dt:tstim(i,2)/dt);
            stimidx(stimidx>length(tspan))=[];
            y(stimidx)=amp(i).*chirp(tspan(stimidx)/1000,0,tstim(i,2)/1000,20);
        end

    case 'pulse'
        y=zeros(size(tspan));% resting state
        for i=1:size(tstim,1)
            % use cumsum to avoid distort bring by round
%             stimidx=round(linspace(tstim(i,1)/dt,(tstim(i,1)+N/intensity*1000)/dt,N));

            stimidx=cumsum(repmat(round(1/intensity*1000/dt),1,N));
            stimidx=stimidx-stimidx(1)+ceil(tstim(i,1)/dt);
            stimidx(stimidx>length(tspan))=[];
            y(stimidx)=1;
        end
        
    case 'poisson'
        y=zeros(size(tspan));% resting state
        for i=1:size(tstim,1)
            dura=tstim(i,2)-tstim(i,1);
            N=round(intensity/1000*dura);
            isis =-log(rand(N,1))/intensity*1000; %inverse transform sampling
            spkt=cumsum (isis);
            if ~isempty(spkt)
                stimidx =floor((spkt-spkt(1)+tstim(i,1))./dt) ;
            else
                stimidx =floor((spkt+tstim(i,1))./dt) ;
            end
            stimidx(stimidx>length(tspan))=[];
            y(stimidx)=1;
        end
end

end
