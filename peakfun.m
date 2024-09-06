
% function for calculating and recording the peak metrics
% arguments = spatial distribution of C, spatial distribution of M, whether
% to return all peaks (summ10=0) or just a summary (summ10=1)
% note can just flip C and M to look at M peaks instead
function [npks0,npks, pklambdas,pkwidths,pkproms,pkheights] = peakfun(Cvals,Mvals,summ10, xset, pkthresh, b1, b2)

[pks,locs,widths,proms] = findpeaks(Cvals, xset, 'MinPeakProminence', pkthresh); % peaks (maxes in C cover)
[pks2,locs2,widths2,proms2] = findpeaks(-Cvals, xset, 'MinPeakProminence', pkthresh); % mins (mins in C cover)
% get the length of any peaks, not saying C needs to be greater than M
% but just the peaks within the region of interest
usepks00 = find(locs <= b2 & locs >= b1);
npks0 = length(pks(usepks00));

usepks = find(Cvals(end, ismember(xset, locs)) > Mvals(end, ismember(xset, locs))); % peaks for which C is greater than M (since these are coral peaks)
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

    %npks = length(usepks); % store total number of peaks. UPDATE: just store the number of peaks in focal region

    % now subset out just the peaks at least 20m from edges
    %usepks = find(locs(usepks) <= b2 & locs(usepks) >= b1); % which peaks are far enough from edges to use
    % intersection of locations in usepks and btw boundaries
    usepks = intersect(find(locs <= b2 & locs >= b1), usepks);
    npks = length(usepks); % store total number of peaks in the focal region

    % and same for the mins
    %usemins = find(locs2(usemins) <= b2 & locs2(usemins) >= b1);
    usemins = intersect(find(locs2 <= b2 & locs2 >= b1), usemins);

    if length(usepks) <= 1 | length(usemins) <= 1% if there was only one peak or fewer
        % calculate the heights, wavelengths, etc.
        %pklambdas = diff(locs2(usemins));
        pklambdas = NaN;
        pkheights = pks(usepks);
        pkwidths = widths(usepks);
        pkproms = proms(usepks);

        %hflag(i) = length(find(abs(pks(usepks)-medh)/medh > etol)); % flag if a peak height is more different than error tolerance
        %lflag(i) = length(find(abs(diff(locs(usepks))-medl)/medl > etol)); % flag if a wavelength is more different than error tolerance

    elseif length(usepks) == 2 % if there only two peaks 
            % calculate the heights, wavelengths, etc.
        pklambdas = diff(locs2(usemins));
        pkheights = pks(usepks);
        pkwidths = widths(usepks);
        pkproms = proms(usepks);

    else % if there are at least 3 peaks
        % remove the edge peaks/mins
        %usepks = usepks(2:(end-1));
        %usemins = usemins(2:(end-1));

        if summ10 ==1
           
        % get the middle 2 peaks (=2 peaks closest to center)
        % pkdist = sort(abs(locs(usepks)));
        % usepks2 = usepks(ismember(abs(locs(usepks)), pkdist(1:2)));
        % usepks2 = usepks2(1:2);% if there are equidistant peaks, both will be included in pkdist(1:2) so need to specify just 2 again here 

        % update: use the middle 3 peaks
        pkdist = sort(abs(locs(usepks)));

        % get the leading mins of middle 3 peaks (=3 peaks closest to center)
        usepks3 = usepks(ismember(abs(locs(usepks)), pkdist(1:3)));
        usepks3 = usepks3(1:3);
        usemins3 = find(ismember(alllocs, locs(usepks3)))-1;
        usemins3 = find(ismember(locs2, alllocs(usemins3))); % get the min indices back in terms of the mins vectors (locs2, pks2, etc.), not alllocs

       % calculate the heights, wavelengths, etc. for the middle 3 peaks
        pklambdas = diff(locs2(usemins3)); % there will only be 2
        pkheights = pks(usepks3);
        pkwidths = widths(usepks3);
        pkproms = proms(usepks3);

        % sort these so 1 = larger wavelength, 2 = smaller wavelength and
        % the other orders are kept consistent
        pklambdas = sort(pklambdas, 'descend'); % sort these on there own since there are only 2
        % update: sort by peak width
        [pkwidths, sortindx] = sort(pkwidths, 'descend'); % update: sort by peak widths
        pkheights = pkheights(sortindx);
        %pkwidths = pkwidths(sortindx);
        pkproms = pkproms(sortindx);

        else

        % calculate all the heights, wavelengths, etc.
        pklambdas = diff(locs2(usemins));
        pkheights = pks(usepks);
        pkwidths = widths(usepks);
        pkproms = proms(usepks);
        end

        %hflag(i) = length(find(abs(pks(usepks)-medh)/medh > etol)); % flag if a peak height is more different than error tolerance
        %lflag(i) = length(find(abs(diff(locs(usepks))-medl)/medl > etol)); % flag if a wavelength is more different than error tolerance

    end

end

