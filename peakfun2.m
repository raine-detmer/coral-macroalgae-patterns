
% function for calculating and recording the max number of peaks (across C
% and M peaks)
% arguments = spatial distribution of C, spatial distribution of M
function [mxpks] = peakfun2(Cvals,Mvals, xset, pkthresh, b1, b2)

% start with coral
[pks,locs,widths,proms] = findpeaks(Cvals, xset, 'MinPeakProminence', pkthresh); % peaks (maxes in C cover)
usepks = find(Cvals(end, ismember(xset, locs)) > Mvals(end, ismember(xset, locs))); % peaks for which C is greater than M (since these are coral peaks)
    % to get C at each peak, need to get the index of xset corresponding to
    % the peak location
    % C values at the final timepoint at the spatial locations where peaks
    % are located: t = end, x = ismember(xset, locs)

    %npks = length(usepks); % store total number of peaks. UPDATE: just store the number of peaks in focal region

    % now subset out just the peaks at least 20m from edges
    %usepks = find(locs(usepks) <= b2 & locs(usepks) >= b1); % which peaks are far enough from edges to use
    % intersection of locations in usepks and btw boundaries
    usepks = intersect(find(locs <= b2 & locs >= b1), usepks);
    npks = length(usepks); % store total number of peaks in the focal region
    npksC = npks;

    % now M
[pks,locs,widths,proms] = findpeaks(Mvals, xset, 'MinPeakProminence', pkthresh); % peaks 
usepks = find(Mvals(end, ismember(xset, locs)) > Cvals(end, ismember(xset, locs))); % peaks for which M is greater than C (since these are coral peaks)
    % to get C at each peak, need to get the index of xset corresponding to
    % the peak location
    % C values at the final timepoint at the spatial locations where peaks
    % are located: t = end, x = ismember(xset, locs)

    %npks = length(usepks); % store total number of peaks. UPDATE: just store the number of peaks in focal region

    % now subset out just the peaks at least 20m from edges
    %usepks = find(locs(usepks) <= b2 & locs(usepks) >= b1); % which peaks are far enough from edges to use
    % intersection of locations in usepks and btw boundaries
    usepks = intersect(find(locs <= b2 & locs >= b1), usepks);
    npks = length(usepks); % store total number of peaks in the focal region
    npksM = npks;



    mxpks = max(npksC, npksM);

end

