
% function for calculating and recording the peak metrics
% arguments = spatial distribution of C (Cvals), spatial distribution of M (Mvals), whether
% to return all peaks (summ10=0) or just the 3 closest to the landscape center (summ10=1), set
% of spatial points (xset) used in simulations that generated Cvals and
% Mvals, min prominence that a peak has to have to count as a patch
% (pkthresh), and lower and upper boundaries of the landscape for peak
% conderation (b1 and b2)
% default is to record metrics of coral peaks, but can just flip order of Cvals and Mvals in the function input to focus on M peaks instead
function [npks0,npks, pklambdas,pkwidths,pkproms,pkheights] = peakfun(Cvals,Mvals,summ10, xset, pkthresh, b1, b2)
% function outputs: npks0 = total number of coral peaks, npks = number of
% coral peaks where coral cover at the peak is greater than macroalgal
% cover, pklambdas = peak wavelengths, pkwidths = peak widths, pkproms =
% peak prominance (height above surrounding coral cover), pkheights = peak
% height (peak height above 0)

[pks,locs,widths,proms] = findpeaks(Cvals, xset, 'MinPeakProminence', pkthresh); % peaks (maxes in C cover)
[pks2,locs2,widths2,proms2] = findpeaks(-Cvals, xset, 'MinPeakProminence', pkthresh); % mins (mins in C cover)
% get the number of C peaks within the region of interest, regardless of whether C is greater than M at
% its peaks
usepks00 = find(locs <= b2 & locs >= b1);
npks0 = length(pks(usepks00));

usepks = find(Cvals(end, ismember(xset, locs)) > Mvals(end, ismember(xset, locs))); % peaks for which C is greater than M (coral patches)
    % to get C at each peak, need to get the index of xset corresponding to
    % the peak location
    % C values at the final timepoint at the spatial locations where peaks
    % are located: t = end, x = ismember(xset, locs)

    % to get leading mins (mins preceding maxes): sort locations of locs(usepks) and locs2, then
    % select the elements of this that are one less than the elements
    % corresponding to the locs(usepks) locations
    alllocs = sort([locs(usepks), locs2]);
    usemins = find(ismember(alllocs, locs(usepks)))-1;
    usemins = usemins(usemins>0); % make sure there were no 0 indices
    usemins = find(ismember(locs2, alllocs(usemins))); % get the min indices back in terms of the mins vectors (locs2, pks2, etc.), not alllocs

    % now subset out just the peaks within the region of interest
    % intersection of locations in usepks and btw boundaries
    usepks = intersect(find(locs <= b2 & locs >= b1), usepks);
    npks = length(usepks); % store total number of peaks in the focal region

    % and do the same for the mins
    usemins = intersect(find(locs2 <= b2 & locs2 >= b1), usemins);

    if length(usepks) <= 1 | length(usemins) <= 1% if there was only one peak or fewer
        % calculate the heights, wavelengths, etc.
        pklambdas = NaN; % no wavelenth since only one peak
        pkheights = pks(usepks);
        pkwidths = widths(usepks);
        pkproms = proms(usepks);

    elseif length(usepks) == 2 % if there only two peaks 
            % calculate the heights, wavelengths, etc.
        pklambdas = diff(locs2(usemins));
        pkheights = pks(usepks);
        pkwidths = widths(usepks);
        pkproms = proms(usepks);

    else % if there are at least 3 peaks

        if summ10 ==1 % if only looking at the 3 middle peaks
           
        
        % order the peak locations
        pkdist = sort(abs(locs(usepks)));

        % get the leading mins of middle 3 peaks (=3 peaks closest to center)
        usepks3 = usepks(ismember(abs(locs(usepks)), pkdist(1:3)));
        usepks3 = usepks3(1:3);
        usemins3 = find(ismember(alllocs, locs(usepks3)))-1;
        usemins3 = find(ismember(locs2, alllocs(usemins3))); % get the min indices back in terms of the mins vectors (locs2, pks2, etc.), not alllocs

       % calculate the heights, wavelengths, etc. for the middle 3 peaks
        pklambdas = diff(locs2(usemins3)); % there will only be 2 wavelength values
        pkheights = pks(usepks3); % heights
        pkwidths = widths(usepks3); % widths
        pkproms = proms(usepks3); % prominances

        % sort for cleaner plotting
        % sort these so 1 = larger wavelength, 2 = smaller wavelength and
        % the other orders are kept consistent
        pklambdas = sort(pklambdas, 'descend'); % sort these on there own since there are only 2
        % sort by peak width
        [pkwidths, sortindx] = sort(pkwidths, 'descend'); % update: sort by peak widths
        pkheights = pkheights(sortindx);
        %pkwidths = pkwidths(sortindx);
        pkproms = pkproms(sortindx);

        else % if not just recording the 3 closes

        % calculate all the heights, wavelengths, etc.
        pklambdas = diff(locs2(usemins));
        pkheights = pks(usepks);
        pkwidths = widths(usepks);
        pkproms = proms(usepks);
        end

      
    end

end

