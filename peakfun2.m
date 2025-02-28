
% function for calculating and recording the max number of peaks (across C
% and M peaks)
% arguments = spatial distribution of C (Cvals), spatial distribution of M (Mvals), set
% of spatial points (xset) used in simulations that generated Cvals and
% Mvals, min prominence that a peak has to have to count as a patch
% (pkthresh), and lower and upper boundaries of the landscape for peak
% conderation (b1 and b2)
function [mxpks] = peakfun2(Cvals,Mvals, xset, pkthresh, b1, b2)

% start with coral
[pks,locs,widths,proms] = findpeaks(Cvals, xset, 'MinPeakProminence', pkthresh); % peaks (maxes in C cover)
usepks = find(Cvals(end, ismember(xset, locs)) > Mvals(end, ismember(xset, locs))); % peaks for which C is greater than M (coral patches)
    % to get C at each peak, need to get the index of xset corresponding to
    % the peak location
    % C values at the final timepoint at the spatial locations where peaks
    % are located: t = end, x = ismember(xset, locs)

    % now subset out just the peaks at are between b1 and b2
    % intersection of locations in usepks and btw boundaries
    usepks = intersect(find(locs <= b2 & locs >= b1), usepks);
    npks = length(usepks); % store total number of peaks in the focal region
    npksC = npks;

    % now repeat fo M
[pks,locs,widths,proms] = findpeaks(Mvals, xset, 'MinPeakProminence', pkthresh); % peaks 
usepks = find(Mvals(end, ismember(xset, locs)) > Cvals(end, ismember(xset, locs))); % peaks for which M is greater than C (macroalgae patches)
    
    % now subset out just the peaks that are between b1 and b2
     % intersection of locations in usepks and btw boundaries
    usepks = intersect(find(locs <= b2 & locs >= b1), usepks);
    npks = length(usepks); % store total number of peaks in the focal region
    npksM = npks;



    mxpks = max(npksC, npksM); % return number of coral peaks and number of macroalgal peaks

end

