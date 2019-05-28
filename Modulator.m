function out = Modulator(binary_data,mode)

%lancia un errore nel caso di input non coerente con il tipo di Modulazione
if mod( length(binary_data),2*mode ) ~= 0
    switch mode
    case 1
        er = 'QPSK';
    case 2
        er = 'QAM16';
    case 3
        er = 'QAM64';
    end
    ModulationError = MException('Modulator:incoerentInput', 'Errore: input non coerente con il tipo di modulazione scelta;\n        Input di lunghezza %d per modulazione di tipo %s.',length(binary_data),er);
    throw( ModulationError )
end

%Tabelle 3GPP TS 36.211 v15.3.0 - Physical channels and modulation
% Paragrafo 7.1
qpsk_iq  = [1+1i; 1-1i; -1+1i; -1-1i]/sqrt(2);
qam16_i = [1 1 3 3 1 1 3 3 -1 -1 -3 -3 -1 -1 -3 -3].'/sqrt(10);
qam16_q = [1 3 1 3 -1 -3 -1 -3 1 3 1 3 -1 -3 -1 -3].'/sqrt(10);
qam64_i = [3 3 1 1 3 3 1 1 5 5 7 7 5 5 7 7 3 3 1 1 3 3 1 1 5 5 7 7 5 5 7 7 -3 -3 -1 -1 -3 -3 -1 -1 -5 -5 -7 -7 -5 -5 -7 -7 -3 -3 -1 -1 -3 -3 -1 -1 -5 -5 -7 -7 -5 -5 -7 -7].'/sqrt(42);  
qam64_q = [3 1 3 1 5 7 5 7 3 1 3 1 5 7 5 7 -3 -1 -3 -1 -5 -7 -5 -7 -3 -1 -3 -1 -5 -7 -5 -7 3 1 3 1 5 7 5 7 3 1 3 1 5 7 5 7 -3 -1 -3 -1 -5 -7 -5 -7 -3 -1 -3 -1 -5 -7 -5 -7].'/sqrt(42);

%Per il tipo di modulazione richiesta viene ristrutturato l'ingresso in una matrice con:
%n colonne pari alla parola della modulazione
%m righe pari al numero di parole da modulare
words_num = length(binary_data)/(2*mode);
in_mat = reshape(binary_data,2*mode,words_num).';      

%Crea una matrice di indici per utilizzare le tabelle
length_array = 2*mode-1 : -1 : 0;
power_array  = 2.^ length_array;
power_array  = repmat(power_array,words_num,1);
index = sum(in_mat.* power_array,2)+1;


switch mode
    case 1
        out = qpsk_iq( index );
    case 2
        out = qam16_i( index ) + 1i*qam16_q( index );
    case 3
        out = qam64_i( index ) + 1i*qam64_q( index );
end
end

