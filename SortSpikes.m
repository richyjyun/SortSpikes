function [spike, snips, p] = SortSpikes(data,fs,stim,range,thresh,ts)
% Summary
%   Returns sorted spikes
%
% Inputs
%   data - filtered data stream or snippets from threshold crossings
%   fs - sampling rate
%   stim - stim times if applicable (for removing noise)
%   thresh - threshold for the snippets (if ndims(data)>1)
%   ts - timestamp for the snippets (if ndims(data)>1)
%
% Outputs
%   spike - array of sorted spike times in seconds
%   snips - matrix of spike traces

if(~exist('stim'))
    stim = [];
end
spikerange = round(-0.0005*fs):1:round(0.0015*fs);

sorttype = 3;
while(sorttype==3)
    close all;
    
    if(any(size(data) == 1))
        % Query threshold
        figure;
        temp = data;
        % if there is stimulation, remove artifact
        if(~isempty(stim))
            stimrange = 0:round(0.001*fs);
            inds = getTrialinds(round(stim*fs),stimrange,length(temp));
            inds(isnan(inds)) = [];
            temp(inds) = nan;
        end
        subplot(1,2,1); plot((1:(30*fs))/fs,temp(1:(30*fs)));
        xlim([0,30]);
        subplot(1,2,2); plot((1:(30*fs))/fs,temp(end-((30*fs)-1):end));
        xlim([0,30]);
        threshold = input('Threshold?');
        close(gcf);
        
        if(threshold >= 0)
            spike = [];
            snips = [];
            return;
        end
        
        % Get snippets
        [snips, cross] = GetSnips(data, threshold, spikerange, 1000, stim,  0.0012, fs);
%         snips = snips-mean(snips);
%         snips = snips-(max(snips)+min(snips))/2;
    else
        snips = data;
        spikerange = range;
        threshold = thresh;
        
        dt = nextEventDT(ts,stim);
        bad = dt<0.0012;
        snips(:,bad) = [];
        ts(bad) = [];
        cross = ts*fs;
    end
    
    % Plot snips
    sorttype = 0;
    while(sorttype==0 || sorttype==3)
        close all;
        figure; plot(snips(:,randsample(size(snips,2),min(size(snips,2),1000))),'k');
        sorttype = input('0) Replot snips, 1) 2 Window Discrim, 2) PCA, 3) Remove noise 4) Re-threshold 5) No spikes');
        if sorttype == 3
            display('Choose noise window')
            line_n = drawline;
            xloc = round(mean(line_n.Position(:,1)));
            
            rm = snips(xloc,:)>min(line_n.Position(:,2)) & snips(xloc,:)<max(line_n.Position(:,2));
            snips(:,rm) = [];
            cross(rm) = [];
        end
        if(sorttype == 5)
            spike = [];
            snips = [];
            return;
        end
    end
    if(sorttype==4), continue; end
end

%% 2 Window discrim
if(sorttype == 1)
    verify = 0;
    while ~verify
        
        display('Choose windows')
        
        avg = nanmean(snips,2);
        
        % Set parameters for 2 window
        p = struct;
        p.thresh = threshold;
        p.threshind = find(spikerange==0);
        
        line_1 = drawline;
        line_1Pos = line_1.Position;
        p.win1del = round(line_1Pos(1))-p.threshind;
        p.win1min = min(line_1Pos(:,2));
        p.win1max = max(line_1Pos(:,2));
        
        line_2 = drawline;
        line_2Pos = line_2.Position;
        p.win2del = round(line_2Pos(1))-p.threshind;
        p.win2min = min(line_2Pos(:,2));
        p.win2max = max(line_2Pos(:,2));
        
        close(gcf)
        
        % Sort
        try
            detected = TwoWindowDiscrim(snips, p);
        catch
            warning('Windows were bad')
            figure; plot(snips(:,randsample(size(snips,2),min(size(snips,2),1000))),'k');
            continue;
        end
        avg = mean(snips(:,detected),2);
        temp = std(snips(:,detected),[],2);
        figure; patch([spikerange/fs,fliplr(spikerange)/fs],...
            [avg+temp;flipud(avg-temp)],'k','facealpha',0.4,'edgealpha',0)
        hold on; plot(spikerange/fs,avg,'k','linewidth',2);
        
        verify = input('Good? (0 if bad, anything else for good)');
        if ~verify
            close(gcf)
            spike = [];
            snips = [];
            return;
        end
    end
    close(gcf);
    
    spike = cross(detected)/fs;
    snips = snips(:,detected);
    
    %% PCA
elseif(sorttype == 2)
    
    % Sort with PCA
    tempsnips = snips;
    tempsnips(isnan(snips)) = 0;
    
    [coeff,score,~,~,explained,~] = pca(tempsnips');
    variance = cumsum(explained/sum(explained));
    r = find(variance > 0.90,1);
    
    Vr = score(:,1:r);
    
    figure; scatter3(Vr(:,1),Vr(:,2),Vr(:,3),'.');
    
    verify = 0;
    while ~verify
        
        in = input('1) Manual boundary, 2) kmeans');
        if in == 1
            
            close(gcf);
            figure; scatter(Vr(:,1),Vr(:,2),'.');
            
            display('Draw boundary')
            e = drawellipse;
            
            in = input('Adjust ellipse and hit enter to continue');
            
            p = CheckEllipse(e.Center(1),e.Center(2),e.SemiAxes(1),e.SemiAxes(2),e.RotationAngle-180, Vr(:,1),Vr(:,2));
            
            detected = p<=1;
            
            avg = mean(snips(:,detected),2);
            temp = std(snips(:,detected),[],2);
            figure; patch([spikerange/fs,fliplr(spikerange)/fs],...
                [avg+temp;flipud(avg-temp)],'k','facealpha',0.4,'edgealpha',0)
            hold on; plot(spikerange/fs,avg,'k','linewidth',2);
            
            verify = input('Good? (0 if bad, anything else for good)');
            if verify
                spike = cross(detected)/fs;
                snips = snips(:,detected);
            else
                close(gcf)
                spike = [];
                snips = [];
                return;
            end
            close(gcf);
            
            %             % For debugging
            %             figure; scatter(Vr(detected,1),Vr(detected,2),'.')
            %             hold on; scatter(Vr(~detected,1),Vr(~detected,2),'r.')
            
        else
            
            ns = input('How many spikes?');
            
            idx = kmeans(Vr,ns,'Replicates',50,'MaxIter',1000);
            figure; subplot(1,2,1);
            for i = 1:max(idx)
                hold on; scatter3(Vr(idx==i,1),Vr(idx==i,2),Vr(idx==i,3),'.');
            end
            subplot(1,2,2); colors = get(gca,'colororder');
            for i = 1:max(idx)
                avg = mean(snips(:,idx==i),2);
                temp = std(snips(:,idx==i),[],2);
                hold on; patch([spikerange/fs,fliplr(spikerange)/fs],...
                    [avg+temp;flipud(avg-temp)],colors(i,:),'facealpha',0.4,'edgealpha',0)
                hold on; plot(avg,'color',colors(i,:),'linewidth',2);
            end
            
            verify = input('Which spike? (0 to change number of spikes)');
            
            if verify
                detected = idx==verify;
                spike = cross(detected)/fs;
                snips = snips(:,detected);
            else
                close(gcf)
                spike = [];
                snips = [];
                return;
            end
        end
        
    end
    close gcf;
end

end