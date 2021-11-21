function [trialinds,bad] = getTrialinds(trig, range, lim)

if(size(trig,1) > size(trig,2))
    trig = trig';
end

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
bad = floor(trialinds(1,:))<=0 | floor(trialinds(end,:))>lim;
trialinds(:,bad) = [];

end