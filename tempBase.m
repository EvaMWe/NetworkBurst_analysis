% calculates the baseline and the temporary calculation of deltaF/F
%
% SYNTAX
% [dff, td_baselind]= tempBase(trace, f, doPlot)
%
% DESCRIPTION
% input
%      trace is the signal trace // vector
%      f is the frequency/frame rate in Hz // scalar
%      doPlot indicates if the result should be displayed // 0 = no, 1 = yes
%
% output
%      dff = the normalized trace according to deltaF/F trace

function [dff,td_baseline] = tempBase( trace, f, doPlot )

nFrames=size(trace,2); %Anzahl Bilder

t1_s=2;
t2_s=5;

t1=round((2*ceil(t1_s/2)+1)*f);
t2=round(t2_s*f);

% Calculate the time dependent baseline
td_baseline=zeros(size(trace));
dff = td_baseline;
F_bar = td_baseline;
F_bar(1,:)=smooth(trace(1,:),t1) ;

for frame=1:nFrames
    if frame <= t2
        [median,~] = Minimum_median(F_bar(1,1:t2-1),5,'Type', 'percent','Dimension',2);
        td_baseline(1,frame) = median;
    else
        [median,~]  = Minimum_median(F_bar(1,frame-t2:frame),5,'Type', 'percent','Dimension',2);
        td_baseline(1,frame) = median;
    end
end
% 
% for frame=1:n_frames
%     rel_f(1,frame)=(meanstack(1,frame)-td_baseline(1,frame))/td_baseline(1,frame);
% end
if doPlot == 1
    figure
    subplot (3,1,1)
    title 'Raw Data';
    plot (meanstack)
    
    subplot(3,1,2)
    title 'Baseline';
    plot (td_baseline)
    
    subplot(3,1,3)
    title 'dff';
    plot(dff)
end
%
end
% 
% 
