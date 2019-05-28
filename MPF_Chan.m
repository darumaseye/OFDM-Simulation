function [zp_received_data,cp_received_data, inv_matrix_zp, inv_matrix_cp] = MPF_Chan(zp_sended_data,cp_sended_data,fading_length,n_of_subcarriers)


%discrete_transfer_func = 0.8*(randn(fading_length,1) + 1i*randn(fading_length,1));

discrete_transfer_func=[0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0 0.21 0.03 0.07].';

channel_matrix = toeplitz([discrete_transfer_func; zeros(length(zp_sended_data)-fading_length,1)],[discrete_transfer_func(1) zeros(1,length(zp_sended_data)-1)]);

zp_received_data = reshape(channel_matrix * zp_sended_data,[],1);
cp_received_data = reshape(channel_matrix * cp_sended_data,[],1);

inv_matrix_cp=inv(diag(fft(discrete_transfer_func,n_of_subcarriers)));
%inv_matrix_zp = inv_matrix_cp;
inv_matrix_zp = pinv(channel_matrix);


end