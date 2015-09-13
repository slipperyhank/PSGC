function X = periodogram(x, N)

window = 'hamming';

if strcmp(window, 'none')
    W = ones(1, N);
elseif strcmp(window, 'hamming')
    W = hamming(N)';
end


K = floor(length(x) / (N / 2)) - 1;

X = zeros(1, N / 2 + 1);

for i=1:K
    y = x(((i - 1) * (N / 2) + 1):((i + 1) * (N / 2))) .* W;
    psdY = fft(y);
    if mod(N, 2) == 0
        psdY = psdY(1:N / 2 + 1);
        psdY = (1 / N) * abs(psdY) .^ 2;
        psdY(2:end - 1) = 2 * psdY(2:end - 1);
        X = X + psdY;
    end
end

X = X / K;
