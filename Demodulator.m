function binary_data = Demodulator(modulated_data,mode)

%array di supporto per la fase decisionale, rappresenta la modulazione
%nel grafico Q,I
qpsk_array  = [ -1 1 ]/sqrt(2);
qam16_array = [ -3 -1 1 3 ]/sqrt(10);
qam64_array = [ -7 -5 -3 -1 1 3 5 7 ]/sqrt(42);

%Tabelle 3GPP TS 36.211 v15.3.0 - Physical channels and modulation
% Paragrafo 7.1
qpsk_i  = [1 1 -1 -1]/sqrt(2);
qpsk_q  = [1 -1 1 -1]/sqrt(2);
qam16_i = [1 1 3 3 1 1 3 3 -1 -1 -3 -3 -1 -1 -3 -3].'/sqrt(10);
qam16_q = [1 3 1 3 -1 -3 -1 -3 1 3 1 3 -1 -3 -1 -3].'/sqrt(10);
qam64_i = [3 3 1 1 3 3 1 1 5 5 7 7 5 5 7 7 3 3 1 1 3 3 1 1 5 5 7 7 5 5 7 7 -3 -3 -1 -1 -3 -3 -1 -1 -5 -5 -7 -7 -5 -5 -7 -7 -3 -3 -1 -1 -3 -3 -1 -1 -5 -5 -7 -7 -5 -5 -7 -7].'/sqrt(42);  
qam64_q = [3 1 3 1 5 7 5 7 3 1 3 1 5 7 5 7 -3 -1 -3 -1 -5 -7 -5 -7 -3 -1 -3 -1 -5 -7 -5 -7 3 1 3 1 5 7 5 7 3 1 3 1 5 7 5 7 -3 -1 -3 -1 -5 -7 -5 -7 -3 -1 -3 -1 -5 -7 -5 -7].'/sqrt(42);

%vettori per l'uso della funzione qamdemod
qpsk_vect = [2 3 0 1];
qam16_vect =[11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
qam64_vect = [47,46,42,43,59,58,62,63,45,44,40,41,57,56,60,61,37,36,32,33,49,48,52,53,39,38,34,35,51,50,54,55,7,6,2,3,19,18,22,23,5,4,0,1,17,16,20,21,13,12,8,9,25,24,28,29,15,14,10,11,27,26,30,31];

switch mode
    case 1
    decisional_array = qpsk_array;
    tab_real = qpsk_i;
    tab_imag = qpsk_q;
    symb_order=qpsk_vect;
    coeff=2;

    case 2
    decisional_array = qam16_array;
    tab_real = qam16_i;
    tab_imag = qam16_q;
    symb_order=qam16_vect;
    coeff=10;
    
    case 3
    decisional_array = qam64_array;
    tab_real = qam64_i;
    tab_imag = qam64_q;
    symb_order=qam64_vect;
    coeff=42;
    
end

if 0
 
i=0;
q=0;
indexes= zeros(1,length(modulated_data));

%Per ogni simbolo rumoroso noisy_sym si calcola parte reale e immaginaria, 
%e per ognuna di esse si trovano le componenti in fase I e in quadratura Q
%associate al simbolo di modulazione.
 k=0;
 for noisy_sym = modulated_data.'
    reale = real(noisy_sym);
    imagin = imag(noisy_sym);
  
    for i = decisional_array

        if reale <= i break;
        end
    end      

    for q = decisional_array

        if imagin <= q break;
        end
    end      
    
  % Usando tab_real/imag come un dizionario, si calcola l'indice associato
  % alle componenti I/Q trovate,
  % tale indice rappresenta la parola di modulazione in decimale 
    
    real_tab_indexes = find(tab_real == i);
    imag_tab_indexes = find(tab_imag == q);
    new_index = intersect(real_tab_indexes,imag_tab_indexes);
    k=k+1;
    indexes(k) = new_index;
 end
 
%Tali indici vengono convertiti in binario e ristrutturati "serialmente" 
columns = 2*mode;
decimal_sym = indexes-1;
binary_data1 = de2bi(decimal_sym.','left-msb'); %2*mode colonne
binary_data = reshape(binary_data1.',1,[]);

else
    
    
modulated_data=modulated_data*sqrt(coeff);   
binary_data = qamdemod(modulated_data,2.^(2*mode),symb_order);
binary_data = de2bi(binary_data,'left-msb');
binary_data = reshape(binary_data.',1,[]);

end


end
