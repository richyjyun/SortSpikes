% Gets snips and time of threshold crossing from data using threshold. 
% Removes any with too large deflections or too close to stim if indicated
function [snips,cross] = GetSnips(data, threshold, range, lim, stim,  dt, fs)

% Find threshold crossing
cross = data < threshold;
cross = find(diff(cross) == 1)';

if(isempty(cross))
    snips = [];
    cross = [];
    return;
end

% Find all snippets
[spiketrials, bad] = getTrialinds(cross,range,length(data));

snips = data(spiketrials);
cross(bad) = [];

if(~isempty(lim))
    % Remove ones with massive deflections
    bad = any(snips<-lim | snips>lim);
    snips(:,bad) = [];
    cross(bad) = [];
end

if(~isempty(stim) && ~isempty(dt) && ~isempty(fs))
    % Remove ones that are too close to stim and clean up noise.
    ind = discretize(cross/fs,stim);
    bad = find(~isnan(ind));
    temp1 = cross(bad)/fs; temp2 = stim(ind(bad));
    lag = temp1(:)-temp2(:);
    bad = bad(lag<dt);
    
    snips(:,bad) = [];
    cross(bad) = [];
end

end