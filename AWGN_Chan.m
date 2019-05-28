%Aggiunge rumore al segnale in ingresso, a seconda del rapporto SNR
%specificato (supposto in dB), se NormSigPow è false calcola la potenza
%del segnale ingresso altrimenti, considera tale potenza normalizzata.
function out = AWGN_Chan(sended_data,NormSigPow,snr)

if NormSigPow
    signal_power = 1;
else
    signal_power = mean(conj(sended_data) .* sended_data);
end
lin_snr =10^(snr/10);
noise_sdp = signal_power / lin_snr;

noise = sqrt(noise_sdp/2)* (randn(length(sended_data),1) + 1i*randn(length(sended_data),1));

out = sended_data + noise;

end