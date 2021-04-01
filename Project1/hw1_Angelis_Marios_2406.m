%Advanced telecommunication systems
%Angelis Marios, AEM:2406
%Implementation of a digital baseband telecommunication system

clear all;
close all;
clc;
pkg load communications;

%We send the same 10000 bits for each SNR. The oversampling factor (OS) is set to 5.
SNR_db=[0,2,4,6,8,10,12,14,16,18,20];
M_vector=[2,4,8,16];
Eb=1;
OS=5;
input_length=10000;
delta_function_tx_filter=sqrt(1/OS)*ones(1,OS);
BER=[];
theoretical_BER=[];

input_initial=randi([0,1],1,input_length);


for d=1:1:length(M_vector)
  M=M_vector(d);

  %Remove the extra bits from the input vector in order the vector to fit exactrly with this modulation 
  input=input_initial(1:log2(M)*fix(input_length/log2(M)));
  input_length=length(input);

  figure(1);
  stem(input);
  title('Input Vector');
  
  %Symbol Encoder
  symbols=[];
  index=1;

  for i=1:log2(M):input_length
    input_bits=input(i:i+log2(M)-1);
    bit_to_int = bin2dec(int2str(input_bits));
    angle=pi/(M/2) + (bit_to_int*pi)/(M/2);  
    symbols(index)=sqrt(log2(M)*Eb)*(exp(j*angle));
    index+=1;
  endfor


  real_values_vector=real(symbols);
  imaginary_values_vector=imag(symbols);

  %Tx filter 
  %Get rid of approximations
  %Clear imaginary part
  for i=1:1:length(imaginary_values_vector)
    if(abs(imaginary_values_vector(i)) <= 3.4837e-16);
      imaginary_values_vector(i)=0;
    endif
  endfor

  %Clear real part
  for i=1:1:length(real_values_vector)
    if(abs(real_values_vector(i)) <= 3.4837e-16);
      real_values_vector(i)=0;
    endif
  endfor

  figure(2);
  subplot(2,1,1);
  stem(real_values_vector);
  title('Real values after symbol encoding');

  subplot(2,1,2);
  stem(imaginary_values_vector);
  title('Imaginary values after symbol encoding');

  real_tx_out=[];
  imaginary_tx_out=[];

  %Perform the convolution for the real and the imaginary part in the tx filter
  index=1;
  for i=1:1:length(symbols)
    real_tx_out(index:index+OS-1)=conv(real_values_vector(i),delta_function_tx_filter);
    imaginary_tx_out(index:index+OS-1)=conv(imaginary_values_vector(i),delta_function_tx_filter);
    index+=OS;
  endfor

  figure(3);
  subplot(3,1,1);
  stem(real_tx_out);
  title('Real values after tx filter');

  subplot(3,1,2);
  stem(imaginary_tx_out);
  title('Imaginary values after tx filter');

  %Iterate through all available SNR values
  for k=1:1:length(SNR_db)
    
    %Compute the new linear SNR value 
    SNR_lin=10^(SNR_db(k)/10);
    noise_variance=Eb/SNR_lin;
    
    %AWGN insertion
    real_awgn_out=real_tx_out + sqrt(noise_variance)*randn(1,length(real_tx_out));
    imaginary_awgn_out=imaginary_tx_out + sqrt(noise_variance)*randn(1,length(imaginary_tx_out));

    figure(4);
    subplot(4,1,1);
    stem(real_awgn_out);
    title('Real values after AWGN insertion');

    subplot(4,1,2);
    stem(imaginary_awgn_out);
    title('Imaginary values after AWGN insertion');

    %Rx filter
    real_rx_out=[];
    imaginary_rx_out=[];

    %Perform the convolution for the real and the imaginary part in the rx filter
    real_rx_out=conv(real_awgn_out,delta_function_tx_filter);
    imaginary_rx_out=conv(imaginary_awgn_out,delta_function_tx_filter);

    figure(5);
    subplot(5,1,1);
    stem(real_rx_out);
    title('Real values after Rx filter');

    subplot(5,1,2);
    stem(imaginary_rx_out);
    title('Imaginary values after Rx filter');

    %Sampling
    %Sampling every (2*OS) discrete points (at the peak of each convolution)
    %Sample for the real and for the imaginary part
    real_part_sample=real_rx_out(OS:OS:length(real_rx_out));
    imaginary_part_sample=imaginary_rx_out(OS:OS:length(imaginary_rx_out));
    
    %Create the total complex decoded symbols vector
    decoded_symbol_vector=real_part_sample + j*imaginary_part_sample;

    figure(6);
    subplot(6,1,1);
    stem(real_part_sample);
    title('Real values after sampling');

    subplot(6,1,2);
    stem(imaginary_part_sample);
    title('Imaginary values after sampling');

    figure(7);
    scatter(real_part_sample,imaginary_part_sample,color='r');
    title('Receiver constellation diagram');

    %Decision device
    decoded_bits=[];
    index=1;
    for i=1:1:length(decoded_symbol_vector)
      
      min_distance=1000000;
      
      %Compare the complex decoded symbol at position i of the decoded_symbol_vector with every possible symbol of this modulation
      for input_bit=0:1:M-1
        angle=pi/(M/2) + (input_bit*pi)/(M/2);
        symbol=sqrt(log2(M)*Eb)*(exp(j*angle));
        
        real_symbol=real(symbol);
        imag_symbol=imag(symbol);
        
        %Get rid of approximations
        %Clear real part
        if(abs(real_symbol) <= 3.4837e-16)
          real_symbol=0;
        endif
        %Clear imaginary part
        if(abs(imag_symbol) <= 3.4837e-16)
          imag_symbol=0;
        endif
        
        %This is one of the symbol of this specific modulation
        symbol=real_symbol + j*imag_symbol;
       
        %Compute the norm between the symbol of this specific modulation and the decoded symbol
        distance=abs(norm(symbol-decoded_symbol_vector(i)));
        
        %Find the symbol with the minimum distance
        if(distance<min_distance);
          decoded_bits(index:index+log2(M)-1)=de2bi(input_bit,log2(M),'left-msb');
          min_distance=distance;
        endif
   
      endfor
      index+=log2(M);
    endfor

    figure(8);
    stem(decoded_bits);
    title('Decoded bits');

    %Compute the actual BER and the theoretical BER for this SNR
    diff=sum(abs(decoded_bits-input));
    BER(d,k)=diff/length(input);
    theoretical_BER(k)=qfunc(sqrt(SNR_lin));
  endfor
endfor

BER;
theoretical_BER;

%Plot the simulation BER for each MPSK modulation in the same graph
figure(9);
semilogy(SNR_db,BER(1,:),'r-*');
hold on;
semilogy(SNR_db,BER(2,:),'b-*');
semilogy(SNR_db,BER(3,:),'g-*');
semilogy(SNR_db,BER(4,:),'c-*');
hold off;
title('Simulation BER for each MPSK modulation');
legend('Simulation BER BPSK', 'Simulation BER QPSK','Simulation BER 8PSK','Simulation BER 16PSK');

%Plot the theoretical BER and the actual BER in the same graph for each MPSK modulation
figure(10);
semilogy(SNR_db,BER(1,:),'r-*');
hold on;
semilogy(SNR_db,theoretical_BER,'b-*');
hold off;
title('Simulation BER vs theoretical BER for BPSK modulation');
legend('Simulation BER BPSK','Theoretical BER BPSK');


