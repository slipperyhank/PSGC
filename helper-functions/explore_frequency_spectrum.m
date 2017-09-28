function [] = explore_frequency_spectrum(data, sampling_rate)

frequency = 1 : (sampling_rate / 2 + 1);
n_channels = size(data, 1);
figure()
for channel = 1:n_channels
    psd = periodogram(data(channel, :), sampling_rate);
    plot(frequency, log(psd), 'linewidth', 2)
    pause()
end
