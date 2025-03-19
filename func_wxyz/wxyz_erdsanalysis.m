function erds = wxyz_erdsanalysis(data, frequency, baseline)
% WXYZ_ERDSANALYSIS This function performs ERD/ERS processing on MEG data
% 
% This function takes data, frequency, and baseline as inputs. Return the
% calculated ERDS as output.
%
% data      - It is used to calculate the raw data of ERD/ERS. The format
%             must be the result obtained after ft_preprocessing. 
% frequency - The ERD/ERS frequency band, in Hz, needs to be calculated.
% baseline  - The baseline segment used in ERD/ERS calculations, in seconds.
% 
% example: 
%   [erds] = wxyz_erdsanalysis(data, frequency, baseline);
% Author: wxyz
% Version: 1.0
% Last revision date : 2024-04-01


% do the general setup of the function
ft_defaults

% Check data
if ~isfield(data, 'time') || ~isfield(data, 'trial')
    error('The data should contains a field name of ''time'' and ''trial''');
end

% Frequency filter
datatmp = data.trial;
method = 'bandpass'; % bandpass / FFT / 
bpopt = [];
bpopt.order = 6;
switch method
    case 'bandpass'
        datatmpfilt = cellfun(@(x) ft_preproc_bandpassfilter(x, data.fsample, frequency, bpopt.order), datatmp, 'UniformOutput', false);
    case 'fft'

end
% Smooth param
smoothWinLength = 200;

nTrial = numel(data.trial);
% Calc trial mean
trialmean = mean(reshape(cell2mat(datatmpfilt), [size(datatmpfilt{1}), nTrial]), 3);

% Calc trial power
datatmppower = cellfun(@(x) (x-trialmean).^2, datatmpfilt, 'UniformOutput', false);

% Calc trial power mean
trialpow = sum(reshape(cell2mat(datatmppower), [size(datatmpfilt{1}), nTrial]), 3)/(nTrial - 1);
clearvars datatmppower

% Smoothing
trialpow = smoothdata(trialpow, 2, 'movmean', smoothWinLength);

% Calc baseline power
bslineidx = [knnsearch(data.time{1}', baseline(1)) knnsearch(data.time{1}', baseline(2))];
bslinepow = mean(trialpow(:, bslineidx), 2);

% Baseline correction 
avgmat = (trialpow - bslinepow) ./ bslinepow * 100;

% collect the results
erds            = keepfields(data, {'label', 'trialinfo', 'fsample', 'grad'});
erds.time       = data.time{1};
erds.avg        = avgmat;
erds.dimord     = 'chan_time';
erds.config.method  = method;
erds.config.param   = baseline;
erds.config.smooth  = smoothWinLength;