%Advanced telecommunication systems
%Angelis Marios, AEM:2406
%Implementation of a digital baseband telecommunication system using OFDM transmission

clear all;
close all;
clc;
graphics_toolkit ("gnuplot")
pkg load communications;

%We send the same 10000 bits for each SNR. The oversampling factor (OS) is set to 1.
SNR_db=[0,2,4,6,8,10,12,14,16,18,20];
M_vector=[16];
Eb=1;
Rb=1000;
OS=1;
input_length=10000;
delta_function_tx_filter=sqrt(1/OS)*ones(1,OS);
ofdm_subcarriers=64;
cyclic_prefix=3;

%Modulation settings
%modulation="MPSK";
modulation="MQAM";

%Channel settings
%channel_type="LTI";
channel_type="AWGN";
h1=0.9+0.9*j;
h2=0.6+0.6*j;
h3=0.3+0.3*j;
h=[h1 h2 h3]

H=fft(h,ofdm_subcarriers);

BER=[];
theoretical_BER=[];

%Bit vector
input_initial=randi([0,1],1,input_length);

for d=1:1:length(M_vector)
  M=M_vector(d);

  %Remove the extra bits from the input vector in order the vector to fit exactrly with this modulation 
  input1=input_initial(1:log2(M)*fix(input_length/log2(M)));
  %length(input1)
  fix((length(input1)/log2(M))/ofdm_subcarriers)
  input=input1(1:log2(M)*ofdm_subcarriers*fix((length(input1)/log2(M))/ofdm_subcarriers));
  input_length=length(input)

  figure(1);
  stem(input);
  title('Input Vector');
  
  %Symbol Encoder
  symbols=[];
  ofdm_out=[];
  equalizer_out=[];
  index=1;

  for i=1:log2(M):input_length
    input_bits=input(i:i+log2(M)-1);
    bit_to_int = bin2dec(int2str(input_bits));
    
    if strcmp(modulation,"MQAM")==1
      real_part = 2 .* floor (bit_to_int ./ (sqrt(M))) - sqrt(M) + 1;
      imaginary_part = -2 .* mod (bit_to_int, (sqrt(M))) + sqrt(M) - 1;
      %Normalized QAM
      %symbols(index) = (Eb/sqrt(10))*(real_part + j.*imaginary_part);
      %Not normalized QAM
      symbols(index) = (real_part + j.*imaginary_part);
      index+=1;
      
    elseif strcmp(modulation,"MPSK")==1
      angle=pi/(M/2) + (bit_to_int*pi)/(M/2);  
      symbols(index)=sqrt(log2(M)*Eb)*(exp(j*angle));
      index+=1;
    endif
  endfor
  
  real_values_vector=real(symbols);
  imaginary_values_vector=imag(symbols);

  %Tx filter 
  %Get rid of approximations
  %Clear imaginary part
  for i=1:1:length(imaginary_values_vector)
    if(abs(imaginary_values_vector(i)) <= 8.4837e-16);
      imaginary_values_vector(i)=0;
    endif
  endfor

  %Clear real part
  for i=1:1:length(real_values_vector)
    if(abs(real_values_vector(i)) <= 8.4837e-16);
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
  for i=1:1:length(real_values_vector)
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
  
  total_tx_out=real_tx_out + j*imaginary_tx_out;
  
  %Iterate through all available SNR values
  for k=1:1:length(SNR_db)
    
    %OFDM encoder
    %OFDM IFFT
    for i=1:ofdm_subcarriers:length(total_tx_out)
      ofdm_idft_out(i:i+ofdm_subcarriers-1)=ifft(total_tx_out(i:i+ofdm_subcarriers-1),ofdm_subcarriers);   
    endfor
  
    ofdm_with_cyclic_prefix=[];
    index=1;
    
    %Add cyclic prefix
    for i=1:ofdm_subcarriers:length(total_tx_out)
      ofdm_with_cyclic_prefix(index:index+ofdm_subcarriers-1)=ofdm_idft_out(i:i+ofdm_subcarriers-1);
      index+=ofdm_subcarriers;
      ofdm_with_cyclic_prefix(index:index+cyclic_prefix-1)=ofdm_idft_out(i:i+cyclic_prefix-1);
      index+=cyclic_prefix;
    endfor

    if strcmp(channel_type,"AWGN")==1
      channel_out=ofdm_with_cyclic_prefix;
    else
      %Flat fading LTI channel
      index=1;
      for i=1:ofdm_subcarriers+cyclic_prefix:length(ofdm_with_cyclic_prefix)
        channel_out(index:index+ofdm_subcarriers+cyclic_prefix+length(h)-2)=conv(ofdm_with_cyclic_prefix(i:i+ofdm_subcarriers+cyclic_prefix-1),h);
        index+=(ofdm_subcarriers+cyclic_prefix+length(h)-1);
      endfor
    endif
    
    %Compute the new linear SNR value 
    SNR_lin=10^(SNR_db(k)/10);
    noise_variance=Eb/SNR_lin;
    
    %AWGN insertion
    awgn_out=awgn(real(channel_out)+imag(channel_out)*j,SNR_db(k),'measured');
    real_awgn_out=real(awgn_out);
    imaginary_awgn_out=imag(awgn_out);

    awgn_out=real_awgn_out + j*imaginary_awgn_out;
    ofdm_dft_out=[];
    
    %OFDM decoder
    awgn_out_without_cyclic_prefix=[];
    
    %Remove cyclic prefix
    if strcmp(channel_type,"AWGN")==1
      index=1;
      for i=1:ofdm_subcarriers+cyclic_prefix:length(awgn_out)
        awgn_out_without_cyclic_prefix(index:index+ofdm_subcarriers-1)=awgn_out(i:i+ofdm_subcarriers-1);
        index+=ofdm_subcarriers;
      endfor
    else
      index=1;
      for i=1:ofdm_subcarriers+cyclic_prefix+length(h)-1:length(awgn_out)
        awgn_out_without_cyclic_prefix(index:index+ofdm_subcarriers-1)=awgn_out(i:i+ofdm_subcarriers-1);
        index+=ofdm_subcarriers;
      endfor
    endif
    
    %OFDM DFT 
    for i=1:ofdm_subcarriers:length(awgn_out_without_cyclic_prefix)
      ofdm_dft_out(i:i+ofdm_subcarriers-1)=fft(awgn_out_without_cyclic_prefix(i:i+ofdm_subcarriers-1),ofdm_subcarriers);
    endfor
    
    %OFDM equalizer
    if strcmp(channel_type,"AWGN")==1
      equalizer_out=ofdm_dft_out;
    else
      for i=1:ofdm_subcarriers:length(ofdm_dft_out)
        equalizer_out(i:i+ofdm_subcarriers-1)=ofdm_dft_out(i:i+ofdm_subcarriers-1)./H;
      endfor
    endif
    
    real_equalizer_out=real(equalizer_out);
    imag_equalizer_out=imag(equalizer_out);
    
    %Get rid of approximations
    %Clear imaginary part
    for i=1:1:length(real_equalizer_out)
      if(abs(real_equalizer_out(i)) <= 8.4837e-16)
        real_equalizer_out(i)=0;
      endif
    endfor

    %Clear real part
    for i=1:1:length(imag_equalizer_out)
      if(abs(imag_equalizer_out(i)) <= 8.4837e-16)
        imag_equalizer_out(i)=0;
      endif
    endfor

    %Rx filter
    real_rx_out=[];
    imaginary_rx_out=[];

    %Perform the convolution for the real and the imaginary part in the rx filter
    real_rx_out=conv(real_equalizer_out,delta_function_tx_filter);
    imaginary_rx_out=conv(imag_equalizer_out,delta_function_tx_filter);
    
    figure(5);
    subplot(5,1,1);
    stem(real_rx_out);
    title('Real values after Rx filter');

    subplot(5,1,2);
    stem(imaginary_rx_out);
    title('Imaginary values after Rx filter');

    %Sampling
    %Sampling every (OS) discrete points (at the peak of each convolution)
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
        if strcmp(modulation,"MQAM")==1
          real_part = 2 .* floor (input_bit ./ (sqrt(M))) - sqrt(M) + 1;
          imaginary_part = -2 .* mod (input_bit, (sqrt(M))) + sqrt(M) - 1;
          %Normalized QAM
          %symbol = (Eb/sqrt(10))*(real_part + j.*imaginary_part);
          %Not normalized QAM
          symbol = (real_part + j.*imaginary_part);
          
        elseif strcmp(modulation,"MPSK")==1
          angle=pi/(M/2) + (input_bit*pi)/(M/2);
          symbol=sqrt(log2(M)*Eb)*(exp(j*angle));
        endif
        
        real_symbol=real(symbol);
        imag_symbol=imag(symbol);

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
    BER
    theoretical_BER(k)=qfunc(sqrt(SNR_lin));
  endfor
endfor

%Plot the theoretical BER and the actual BER in the same graph for each MPSK modulation
figure(10);
semilogy(SNR_db,theoretical_BER,'r','linewidth',2)
hold on;
semilogy(SNR_db,BER(1,:),'b','linewidth',2)
hold off;
legend('AWGN Theoretical BER','Simulation BER');