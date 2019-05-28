%% Simulazione di OFDM
% Lo script è strutturato in 3 parti:
%Nella fase 0, viene fatta una semplice verifica dei concetti che si
%trovano alla base della OFDM: 
%  viene mostrata la funzione di pulse shapin spiegando lo scopo di ogni sua parte;
%  viene poi disegnata la funzione di autocorrelazione Pm spiegando come
%   la sua particolare forma è una soluzione per la IBI;
%  processo simile viene fatto per la soluzione della ICI disegnando la
%  correlazione incrociata Rpmpj
%
%Nella fase 1, viene impostata una simulazione costituita da Trasmettitore
%e Ricevitore considerando un semplice canale AWGN (implementato nella
%funione AWGN_Chan) evidenziando il rapporto con la DFT e IDFT. 
%Per mostrare il comportamento da un punto di vista di errore, viene
%calcolato il Bit Error Rate, per fare ciò si fa uso di un Modulatore QAM.
%
%Nella fase 2, viene ripetuta la simulazione della fase 1 utilizzando un
%semplice canale di multipath fading del tipo: somma{Ak * delta(t-Tk)} con k=0,1,...,L
%Vengono quindi spieagte le tecniche di Zero Padding e di Cyclic Prefix
%utilizzando un semplice equalizzatore di tipo Zero Forcing.

%% Inizializzazione dello script (funzione 'clear' e aggiunta della cartella di lancio nel percorso di esecuzione di Matlab)

close all; clear all; clc;
addpath(fileparts(which('OFDM_Simulation.m'))); % Comando per aggiungere la cartella al percorso di ricerca di matlab

%% Fase 0, Verifica di assenza di ISI per l'OFDM

n_of_subcarriers_test = 7;
n_of_blocks_test = 1;
symbol_duration = 66600/999 * 10^-3; 
MTs = n_of_subcarriers_test * symbol_duration;

t = 0:0.01:10;
t2 = -1:0.001:1;
m = 0:1:(n_of_subcarriers_test -1 );

% Pm(t): La funzione sagomatrice dell'OFDM e' composta da 3 parti:
% 1. coefficiente pari al reciproco della norma e serve per normalizzare la
% funzione pm per risolvere la IBI:
norm_coeff = 1 / sqrt(MTs);
% 2. Funzione ortogonale esponenziale per risolvere la ICI nello stesso
% blocco (l = k) e per implementare l-OFDM con la FFT:
orthogonal_func = exp(1i * 2 * pi * m.' * t / MTs ); 
% 3. Funzione di finestramento (ideale) per imporre alla Pm una durata
% inferiore di 2MTs per risolvere la ICI tra blocchi diversi (l ~= k),
% (questa funzione agisce insieme alla normalizzazione anche sulla IBI):
rect_func = ones(1,numel(t)).*( abs( t-MTs ) < MTs/2 );
 
p_m = norm_coeff .* orthogonal_func .* rect_func;
figure('Name','Funzione Sagomatrice della OFDM: Pm'); 
subplot(311); plot(abs(p_m));
title('Grafico della Pm (modulo)')

subplot(312); plot(angle(p_m));
title('Grafico della Pm (fase)')

subplot(313); plot(abs(fft(p_m)));
title(sprintf('Grafico dello spettro della Pm (su %d sottoportanti)',n_of_subcarriers_test));


% Dimostrazione di assenza di IBI.
% Lo scopo e' quello di far vedere che grazie alla funzione di
% finestramento e alla normalizzazione di Pm , Rpmpm((l-k)MTs) = 1 per l=k 
% e Rpmpm((l-k)MTs) = 0 per l~=k e per ogni m.
corr_figure = figure('Name','Autocorrelazione di Pm');
clf;
for i = m
    %sum(p_m .* conj(p_m),2)
    y = xcorr(p_m(i+1,:),'coeff');
    corr_figure;
    hold on;
    plot(t2,y);
end
title('Grafico delle autocorrelazioni di Pm per ogni m: Indice di presenza di IBI');
annotation('textbox', [.7 .26 .2 .24], 'String',"Rpmpm((l-k)MTs)=1 per k=l e Rpmpm((l-k)MTs)=0 per k~=l");

% Dimostrazione di assenza di ICI.
% Lo scopo e' quello di mostrare che la correlazione incrociata della
% funzione Pm e' nulla nei punti di campionamento: Rpmpj((l-k)MTs) = 0.
% Questo grazie alla funzione di finestramento che delimita Rpmpj ai punti
% t= (l-k)MTs con l ~= k (caso ICI da blocchi diversi) e grazie
% all'ortogonalita' della funzione esponenziale che porta a zero il punto
% t= (l-k)MTs con l=k (caso ICI da blocchi uguali).
corr_figure = figure ('Name','Correlazione Incrociata  PmPj');
for i = m
    for j = m
        if i ~= j
            y = xcorr(p_m(i+1,:),p_m(j+1,:),'coeff');
            corr_figure;
            hold on;
            plot(t2,y)
        end
 
    end
end
title('Grafico delle autocorrelazioni incrociate di Pm per ogni m: Indice di presenza di ICI')
annotation('textbox', [.7 .26 .2 .24], 'String',"Rpmpj((l-k)MTs)=0 per k=l e Rpmpj((l-k)MTs)=0 per k~=l");




%% Fase 1, Studio in canale AWGN: INIZIALIZZAZIONE

n_of_subcarriers = 76;
modulation_bits = 6;
mode=modulation_bits/2;
n_of_blocks = 3;
N = n_of_subcarriers;

t = 0:1:N-1;
m = 0:1:(n_of_subcarriers-1);

ber_analysis_mode =true;
normalized_signal_power = true;

if ber_analysis_mode
    snr_vector = 0:50;
    repetition_vector = 1:1:100;
else
    snr_vector = 60;
    repetition_vector = 1;
end

%Preallocation
received_symbols = zeros(n_of_subcarriers * n_of_blocks,1);
sended_data = zeros(n_of_subcarriers * n_of_blocks,1);
mean_vector = zeros(1,length(repetition_vector));

ber_vector = zeros(1,length(snr_vector));


disp('Fase 1a di "AWGN: INIZIALIZZAZIONE", ESEGUITA');

%% GENERAZIONE DI DATI & MODULAZIONE (F1, Tx)

n_of_bits = n_of_subcarriers * modulation_bits * n_of_blocks;
binary_data = randi([0,1],n_of_bits,1).';
modulated_data = Modulator(binary_data,mode);
%dec=bi2de(reshape(binary_data,[modulation_bits,n_of_subcarriers*n_of_blocks]).','left-msb');
%vect=[11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
%modulated_data=qammod(dec,2.^modulation_bits,vect)/sqrt(10);

blocks_of_data = reshape(modulated_data,[n_of_subcarriers, n_of_blocks]);

disp('Fase 1b di "GENERAZIONE DI DATI & MODULAZIONE (F1, Tx)", ESEGUITA');

%% MULTIPLAZIONE OFDM (F1,Tx)
% Definizione matriciale per l'OFDM:
% Il calcolo della DFT in forma matriciale e' dato dal prodotto hermitiano
% del vettore (colonna) x (sequenza nel tempo) e il vettore della base
% ortogonale in e. Da questo si dimostra che la DFT è il prodotto
% matriciale tra la matrice della base ortogonale in e (hermitiana)
% ed il vettore x.

F = exp(1i * 2 * pi * t.' * m / N);
Fh = transpose( conj(F) );
DFT = Fh * modulated_data(1:n_of_subcarriers);

% Considerando F matrice che descrive una base ortogonale per lo spazio dei
% segnali discreti (composta da vettori Fk), vado a rappresentare un 
% generico segnale discreto come combinazione lineare su F.
% Se un coefficiente è dato dalla formula di Fourier: (prodotto scalare 
% tra x e Fk)/(norma di Fk) allora tale coefficiente è uguale a
% DFT(x)/norma quindi:

IDFT = (F * DFT)/N;

%Per verificare, modulated_data == IDFT ... ma qui Matlab impazzisce.

%Per concludere, invece che usare la funzione ifft:
%sended_data = ifft(modulated_data,n_of_subcarriers);
%implemento l'OFDM con una forma matriciale ripetendo l-operazione per ogni blocco:

for i = 1:n_of_blocks
    %sended_data = (F * modulated_data)/n_of_subcarriers;
    sended_data(1+(i-1)*n_of_subcarriers:i*n_of_subcarriers) = (F * blocks_of_data(:,i))/N;

end

disp('Fase 1c di "MULTIPLAZIONE OFDM (F1,Tx)", ESEGUITA');

%% PLOTTING (F1,Tx)
f_trasmission = figure('Name','Transmission in AWGN Channel');
clf;
subplot(311);
stairs(binary_data);
title(sprintf('Bit da Inviare',n_of_subcarriers));

subplot(312);
stem(abs(sended_data));
title(sprintf('Dati Multiplati su %d sottoportanti: IFFT',n_of_subcarriers));

subplot(313);
plot(abs(sended_data));
title('Dati dopo il convertitore D/A: Segnale Inviato');

%scatterplot(modulated_data);
%title('Dati Modulati: Simboli QAM');

disp('Fase 1d di "PLOTTING (F1,Tx)", ESEGUITA');

%% PASSAGGIO NEL CANALE AWGN, DEMULTIPLAZIONE, DEMODULAZIONE (F1,Rx)

k=0;
for snr = snr_vector
 disp(sprintf('Fase 1e Ciclo %d di %d',k+1,length(snr_vector)));
 for i = repetition_vector
    
  received_data = AWGN_Chan(sended_data,normalized_signal_power,snr);
     
  if ber_analysis_mode
           
      blocks_of_data = reshape(received_data,[n_of_subcarriers, n_of_blocks ]);
      received_symbols=fft(blocks_of_data,N);
      
  else %OFDM con matrice,meno efficiente
      received_data = AWGN_Chan(sended_data,0,snr); 
      blocks_of_data = reshape(received_data,[n_of_subcarriers,n_of_blocks]);
      received_symbols = Fh * (blocks_of_data);

  end
  
  serial_symbols = reshape(received_symbols,[],1);
  received_bits = Demodulator(serial_symbols,mode); % Rx De-Modulation
  mean_vector(i) = mean(bitxor(received_bits,binary_data)); % Calcolo del BER


 end
 
 k=k+1;
 ber_vector(k) = mean(mean_vector);

 mean_vector = zeros(1,length(repetition_vector));

end

if ber_analysis_mode  
    figure('Name','BER')
    semilogy(snr_vector,ber_vector,'Marker','s','LineStyle',':','Color','b');
    title(sprintf('BER-Analysis, con SNR da 0 (rumore infinito) a %d',snr_vector(numel(snr_vector))));
    xlabel('SNR (dB)');
    ylabel('BER');
    grid on;

end

disp('Fase 1e di "PASSAGGIO NEL CANALE AWGN, DEMULTIPLAZIONE, DEMODULAZIONE (F1,Rx)", ESEGUITA');

%% PLOTTING (F1, Rx)
    f_reception = figure('Name','Reception  in AWGN Channel');
    clf;
    subplot(311);
    plot(abs(received_data));
    title('Segnale dopo il canale AWGN');

    subplot(312);
    stem(abs(received_data));
    title('Segnale dopo il convertitore A/D: Simboli OFDM');

    subplot(313);
    stairs(received_bits);
    title('Bit ricevuti');
    
    scatterplot(reshape(received_symbols,[],1));
    title('Dati dopo la fft: Simboli QAM rumorosi');

disp('Fase 1f di "PLOTTING (F1, Rx)", ESEGUITA');





%% FASE 2, Studio in un canale awgn con multipath fading: INIZIALIZZAZIONE
n_of_subcarriers = 76;
modulation_bits = 4;
mode=modulation_bits/2;
n_of_blocks = 3;
N = n_of_subcarriers;
fading_length = 11;

t = 0:1:N-1;
m = 0:1:(n_of_subcarriers-1);
F = exp(1i * 2 * pi * t.' * m / N);
Fh = transpose( conj(F) );

ber_analysis_mode =true;
normalized_signal_power = true;

if ber_analysis_mode
    snr_vector = 0:60;
    repetition_vector = 1:50;
else
    snr_vector = 60;
    repetition_vector = 1;
end

%Preallocation
sended_data = zeros(n_of_subcarriers * n_of_blocks,1);

zp_mean_vector = zeros(1,length(repetition_vector));
cp_mean_vector = zeros(1,length(repetition_vector));
zp_ber_vector = zeros(1,length(snr_vector));
cp_ber_vector = zeros(1,length(snr_vector));

disp('Fase 2a di "AWGN+MPF: INIZIALIZZAZIONE", ESEGUITA');

%% GENERAZIONE DI DATI & MODULAZIONE (F2, Tx)

n_of_bits = n_of_subcarriers * modulation_bits * n_of_blocks;
binary_data = randi([0,1],n_of_bits,1).';
modulated_data = Modulator(binary_data,mode);

blocks_of_modulated_data = reshape(modulated_data,[n_of_subcarriers, n_of_blocks]);

disp('Fase 2b di "GENERAZIONE DI DATI & MODULAZIONE (F2, Tx)", ESEGUITA');

%% MULTIPLAZIONE ZP-OFDM & CP-OFDM (F2,Tx)

ofdm_data = (F * blocks_of_modulated_data)/n_of_subcarriers;

zp_matrix = [eye(n_of_subcarriers); zeros(fading_length,n_of_subcarriers)];
cp_matrix = [[zeros(fading_length, n_of_subcarriers - fading_length) eye(fading_length)] ; eye(n_of_subcarriers)];

blocks_of_zp_ofdm_data = zp_matrix * ofdm_data;
blocks_of_cp_ofdm_data = cp_matrix * ofdm_data;

zp_sended_data = reshape(blocks_of_zp_ofdm_data,[],1);
cp_sended_data = reshape(blocks_of_cp_ofdm_data,[],1); 

disp('Fase 2c di "MULTIPLAZIONE OFDM (F2,Tx)", ESEGUITA');

%% PLOTTING (F2,Tx)

figure('Name','Transmission in AWGN+MPF Channel with ZP-OFDM');
clf;
subplot(411);
stairs(binary_data);
title(sprintf('Bit da Inviare',n_of_subcarriers));
subplot(412)
stem(abs(reshape(ofdm_data,[],1)));
title(sprintf('Dati Modulati e Multiplati su %d sottoportanti: IFFT',n_of_subcarriers));
subplot(413);
stem(abs(zp_sended_data));
title('Dati Paddati con ZP');
subplot(414);
stem(abs(cp_sended_data));
title('Dati Paddati con CP');


disp('Fase 2d di "PLOTTING (F2,Tx)", ESEGUITA');

%% PASSAGGIO NEL CANALE AWGN CON M.P.F., EQUALIZZAZIONE, DEMULTIPLAZIONE & DEMODULAZIONE (F2,Rx) 

%Variabili usate per togliere lo zero padding
l=1:fading_length;
Fl = exp(1i * 2 * pi * l.' * m / N);
B= [F; Fl];
Bh = transpose( conj(B) );

%variabili usate per il cyclic prefix
C = [zeros(n_of_subcarriers,fading_length) eye(n_of_subcarriers) ];

%Passaggio nel canale

k=0;
for snr = snr_vector
 disp(sprintf('Fase 2e Ciclo %d di %d',k+1,length(snr_vector)));
 for i = repetition_vector
    
     [zp_faded_data,cp_faded_data, inv_channel_matrix_zp,inv_channel_matrix_cp]= MPF_Chan(blocks_of_zp_ofdm_data,blocks_of_cp_ofdm_data,fading_length,n_of_subcarriers);

     
    zp_received_data = AWGN_Chan(zp_faded_data,normalized_signal_power,snr);
    cp_received_data = AWGN_Chan(cp_faded_data,normalized_signal_power,snr);
    
    zp_received_data = reshape(zp_received_data,n_of_subcarriers+fading_length, n_of_blocks);
    cp_received_data = reshape(cp_received_data,n_of_subcarriers+fading_length, n_of_blocks);
    
    zp_equalized_data = inv_channel_matrix_zp * zp_received_data;
    %cp_equalized_data = inv_channel_matrix_cp * cp_received_data;
     %zp_equalized_data =  zp_received_data;
     cp_equalized_data =  cp_received_data;

    zp_received_symbols = Bh * zp_equalized_data;
    cp_received_symbols = C * cp_equalized_data;
    cp_received_symbols = fft(cp_received_symbols,N);
    
    %equalizz con diag
    %zp_received_symbols = inv_channel_matrix_zp * zp_received_symbols;
    cp_received_symbols = inv_channel_matrix_cp * cp_received_symbols;

    
    serial_zp_received_symbols = reshape(zp_received_symbols,[],1);
    serial_cp_received_symbols = reshape(cp_received_symbols,[],1);

    zp_received_bits = Demodulator(serial_zp_received_symbols,mode); % Rx De-Modulation
    cp_received_bits = Demodulator(serial_cp_received_symbols,mode); % Rx De-Modulation
    
    zp_mean_vector(i) = mean(bitxor(zp_received_bits,binary_data)); % Calcolo del BER
    cp_mean_vector(i) = mean(bitxor(cp_received_bits,binary_data)); 

 end
 
 k=k+1;
 zp_ber_vector(k) = mean(zp_mean_vector);
 cp_ber_vector(k) = mean(cp_mean_vector);
 zp_mean_vector = zeros(1,length(repetition_vector));
 cp_mean_vector = zeros(1,length(repetition_vector));

 
end

 if ber_analysis_mode    
     
    figure('Name','BER for ZP-OFDM & CP-OFDM')
    semilogy(snr_vector,zp_ber_vector,'Marker','s','LineStyle',':','Color','b');
    title(sprintf('BER-Analysis for ZP-OFDM & CP-OFDM, con SNR da 0 (rumore infinito) a %d',snr_vector(numel(snr_vector))));
    hold on;
    semilogy(snr_vector,cp_ber_vector,'Marker','*','LineStyle',':','Color','r');
    legend('zp-ofdm','cp-ofdm');
    xlabel('SNR (dB)');
    ylabel('BER');
    grid on;
 
 end

disp('Fase 2e di "PASSAGGIO NEL CANALE AWGN+MPF, EQUALIZZAZIONE, DEMULTIPLAZIONE & DEMODULAZIONE (F2,Rx)", ESEGUITA');

%% PLOTTING (F2,Rx)
    p=reshape((0:n_of_subcarriers*n_of_blocks-1).',n_of_subcarriers,n_of_blocks);
    p1=reshape((0:(n_of_subcarriers+fading_length)*n_of_blocks-1).',n_of_subcarriers+fading_length,n_of_blocks);

    zp_reception = figure('Name','Reception, ZP-OFDM');
    clf;
    subplot(311);
    plot(p1,abs(zp_received_data));
    title('Segnale dopo il canale con Multipath Fading e AWGN');
    
    subplot(312);
    stem(p1,abs(zp_equalized_data));
    title('Dati Equalizzati');
    
    subplot(313);
    stairs(zp_received_bits);
    title('Bit Ricevuti');
    
    scatterplot(serial_zp_received_symbols);
    title('Flusso di Dati con ZP dopo la fft: Simboli QAM rumorosi');
    
    
    figure('Name','Reception, CP-OFDM');
    clf;
    subplot(311);
    plot(p1,abs(cp_received_data));
    title('Segnale dopo il canale con Multipath Fading e AWGN');

    subplot(312);
    stem(p,abs(cp_received_symbols));
    title('Dati Equalizzati');
    
    subplot(313);
    stairs(cp_received_bits);
    title('Bit Ricevuti');
    
    scatterplot(serial_cp_received_symbols);
    title('Simboli QAM rumorosi dal flusso con Cyclic Prefix');
    
    disp('Fase 2f di "PLOTTING (F2, Rx)", ESEGUITA');
