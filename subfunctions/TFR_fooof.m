function [out] = TFR_fooof(TFR)

clear fspctrm
clear fooofedtrl

for t = 1:length(TFR.time)
    cfg = [];
    cfg.latency = TFR.time(t);
    tmp = ft_selectdata(cfg, TFR);
    
    % Frequency vector (assumes evenly spaced)
    freqs = tmp.freq(:);                 % column
    df = median(diff(freqs));            % frequency resolution (Hz)
    
    % Choose peak width limits: lower bound ~ 2 * df (per FOOOF recommendation)
    lb = 2 * df;
    ub = 12;                             % FOOOF default upper bound
    if ub <= lb
        ub = lb + 1;                     % ensure upper > lower
    end
    
    % FOOOF settings
    settings = struct();
    settings.peak_width_limits = [lb, ub];
    
    f_range = [freqs(1), freqs(end)];
    
    % Preallocate per-time-slice array
    powspctrmff = nan(numel(tmp.label), numel(freqs));
    
    for chan = 1:numel(tmp.label)
        psd = tmp.powspctrm(chan, :).';  % column
        
        % Run FOOOF
        fooof_results = fooof(freqs, psd, f_range, settings, true);
        
        % Reconstructed peaks = full fit minus aperiodic fit
        powspctrmff(chan, :) = fooof_results.fooofed_spectrum - fooof_results.ap_fit;
    end
    
    fspctrm(:, :, t) = powspctrmff;
end
fooofedtrl = fspctrm;


out=TFR;
out.powspctrm = fooofedtrl;