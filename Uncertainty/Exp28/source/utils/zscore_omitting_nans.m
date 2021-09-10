function [zscored_data] = zscore_omitting_nans(x)
    zscored_data = (x - nanmean(x))/nanstd(x);
end