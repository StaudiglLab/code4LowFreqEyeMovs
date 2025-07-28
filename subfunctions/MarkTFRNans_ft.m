function out = MarkTFRNans_ft(data,spike_info)
%% replace the marked artifacts as NaNs
out = data;
for ichan = 1:length(data.label)
    spike_info_chan = spike_info{ichan};
    if ~isempty(spike_info_chan)
        for itrial = 1:size(spike_info_chan,1)
            start_idx = find(spike_info_chan(itrial,2)<=data.time,1,'first');
            end_idx = find(spike_info_chan(itrial,3)>=data.time,1,'last');
            out.powspctrm(spike_info_chan(itrial),ichan,:,start_idx:end_idx)=nan;
        end
    end
end