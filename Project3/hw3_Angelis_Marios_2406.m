%Advanced telecommunication systems
%Angelis Marios, AEM:2406
%Implementation of a digital baseband telecommunication system with a Flat fading LTI channel, Maximum Ratio Combining (MRC) and MIMO with least squares channel inversion (LS)
%Available modulations: MPSK, MQAM

clear all;
close all;
clc;
pkg load communications;
graphics_toolkit ("gnuplot")

%We send the same 10000 bits for each SNR. The oversampling factor (OS) is set to 5.
SNR_db=[0,2,4,6,8,10,12,14,16,18,20];
M_vector=[16];
Eb=1;
Rb=1000;
OS=1;
T=1;
input_length=10000;

%Modulation settings
%modulation="MPSK";
modulation="MQAM";

%MIMO settings
MIMO_on=0;
MIMO_antennas=2;

%MRC settings
diversity_branches=2;
%Equalizer control bit
equalizer_on=1;
%Enable h+c error in the receiver
h_error_on=0;
%Enable approximation of the h in the receiver
h_approximate_on=0;
%Set c variable
c=0;
minimum_transmitter_symbols=0;
delta_function_tx_filter=sqrt(1/OS)*ones(1,OS);
BER=[];
goodput=[];
theoretical_awgn_only_BER=[];
theoretical_flat_fading_BER=[];

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
  for i=1:1:length(real_values_vector)
    %Clear real part
    if(abs(real_values_vector(i)) <= 3.4837e-16)
      real_values_vector(i)=0;
    endif
    %Clear imaginary part
    if(abs(imaginary_values_vector(i)) <= 3.4837e-16)
      imaginary_values_vector(i)=0;
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
    
    %Flat fading channel modeled by one TAP
    h=[];
    total_tx_out=real_tx_out + j*imaginary_tx_out;
    channel_out=[];
    
    %Compute the number of transmitted symbols before the variable h changes value
    number_of_transmitted_symbols=int64(T*Rb/log2(M))
    
    if number_of_transmitted_symbols*OS>=length(total_tx_out)
      number_of_transmitted_symbols=length(total_tx_out)/OS
      minimum_transmitter_symbols=1
    endif
    
    %MIMO
    if MIMO_on==1

      for antenna=1:1:MIMO_antennas
        temp_channel_out=[];
        temp_h=[];
        h_index=0;
        index=1;

        for i=1:1:length(total_tx_out)
          if mod(i,((number_of_transmitted_symbols*OS)))==1
            
            %Fill one row of the matrix H
            for w=1:1:MIMO_antennas
              h_index++;
              temp_h(h_index)=sqrt(1/(sqrt(2)))*(randn(1,1) +j*randn(1,1));
            endfor
          endif
          if mod(i,MIMO_antennas)!=1
            continue
          endif
          
          temp_channel_out(index)=0;
          for w=1:1:MIMO_antennas
            temp_channel_out(index)+= (temp_h(h_index-w+1)*total_tx_out(i+MIMO_antennas-w));
          endfor
          
          index++;
          
        endfor
        h=[h;temp_h];
        channel_out=[channel_out;temp_channel_out];
      endfor
    else
      for index=1:1:diversity_branches
        temp_channel_out=[];
        temp_h=[];
        h_index=0;
        for i=1:1:length(total_tx_out)
          if mod(i,(number_of_transmitted_symbols*OS))==1
            h_index++;
            temp_h(h_index)=sqrt(1/(sqrt(2)))*(randn(1,1) +j*randn(1,1));
          endif
          temp_channel_out(i)=temp_h(h_index)*total_tx_out(i);
        endfor
        h=[h;temp_h];
        channel_out=[channel_out;temp_channel_out];
      endfor
    endif
    
    real_channel_out=real(channel_out);
    imaginary_channel_out=imag(channel_out);
    
    %Compute the new linear SNR value 
    SNR_lin=10^(SNR_db(k)/10);
    noise_variance=Eb/SNR_lin;
    
    %AWGN insertion
    real_awgn_out=real_channel_out + sqrt(noise_variance)*randn(1,length(real_channel_out));
    imaginary_awgn_out=imaginary_channel_out + sqrt(noise_variance)*randn(1,length(imaginary_channel_out));

    total_awgn_out=real_awgn_out +j*imaginary_awgn_out;
    
    %Equalizer
    equalizer_out=[];
    
    if equalizer_on==1
      
      if MIMO_on==1
        h_column=1;
        H=[];
        index=1;
          
        %Perform equalization with multiple diversity branches
        for i=1:1:length(total_awgn_out)
          if minimum_transmitter_symbols==0 && mod(i,((number_of_transmitted_symbols)*OS)/MIMO_antennas)==1
            
            for row=1:1:MIMO_antennas
              for w=1:1:MIMO_antennas
                H(row,w)=h(row,h_column+w-1);
              endfor
            endfor

            h_column+=MIMO_antennas;
            H_inverse=(inv((H')*H))*H';
          elseif minimum_transmitter_symbols==1 && i==1
            
            for row=1:1:MIMO_antennas
              for w=1:1:MIMO_antennas
                H(row,w)=h(row,h_column+w-1);
              endfor
            endfor

            h_column+=MIMO_antennas;
            H_inverse=(inv((H')*H))*H';
          endif
          
          equalizer_out(index:index+MIMO_antennas-1)=H_inverse*total_awgn_out(:,i);
          index+=MIMO_antennas;
        endfor
      else
        %Simple equalization with one diversity branch
        if diversity_branches==1
          h_index=0;
          
          %c error mode is enabled
          if h_error_on==1
            h=h+c;
          endif
          
          %Get rid of the c variable if the h_approximate mode is on
          if h_approximate_on==1
            h(1,:) -= mean(real(h(1,:)));
          endif
            
          %Perform equalization with one diversity branch
          for i=1:1:length(total_awgn_out)
            if mod(i,(number_of_transmitted_symbols*OS))==1
              h_index++;
            endif
            equalizer_out(i)=(total_awgn_out(i)*conj(h(h_index)))/(real(h(h_index))^2 + imag(h(h_index))^2);
          endfor
          
        else
          %Maximum ratio combining (MRC) equalization
          h_index=0;
          
          %Perform equalization with multiple diversity branches
          for i=1:1:length(total_awgn_out)
            if mod(i,(number_of_transmitted_symbols*OS))==1
              h_index++;
              h_appr_norms = 0;
              for ii=1:1:diversity_branches
                h_appr_norms += norm(h(ii,h_index)).^2;
              endfor
              h_inverse=(h(:,h_index)')/h_appr_norms;
            endif
            
            equalizer_out(i)=(h_inverse*total_awgn_out(:,i));
          endfor
        endif
      endif
      
      real_equalizer_out=real(equalizer_out);
      imag_equalizer_out=imag(equalizer_out);
      
      %Get rid of approximations
      for i=1:1:length(equalizer_out)
        %Clear real part
        if(abs(real_equalizer_out(i)) <= 3.4837e-16)
          real_equalizer_out(i)=0;
        endif
        %Clear imaginary part
        if(abs(imag_equalizer_out(i)) <= 3.4837e-16)
          imag_equalizer_out(i)=0;
        endif
      endfor
  
      equalizer_out=real_equalizer_out + j*imag_equalizer_out;
    else
      equalizer_out=real_awgn_out +j*imaginary_awgn_out;
    endif

    %Rx filter
    real_rx_out=[];
    imaginary_rx_out=[];

    %Perform the convolution for the real and the imaginary part in the rx filter
    real_rx_out=conv(real(equalizer_out),delta_function_tx_filter);
    imaginary_rx_out=conv(imag(equalizer_out),delta_function_tx_filter);

    figure(4);
    subplot(4,1,1);
    stem(real_rx_out);
    title('Real values after Rx filter');

    subplot(4,1,2);
    stem(imaginary_rx_out);
    title('Imaginary values after Rx filter');

    %Sampling
    %Sampling every (OS) discrete points (at the peak of each convolution)
    %Sample for the real and for the imaginary part
    real_part_sample=real_rx_out(OS:OS:length(real_rx_out));
    imaginary_part_sample=imaginary_rx_out(OS:OS:length(imaginary_rx_out));
    
    %Create the total complex decoded symbol vector
    decoded_symbol_vector=real_part_sample + j*imaginary_part_sample;

    figure(5);
    subplot(5,1,1);
    stem(real_part_sample);
    title('Real values after sampling');

    subplot(5,1,2);
    stem(imaginary_part_sample);
    title('Imaginary values after sampling');

    figure(6);
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

    figure(7);
    stem(decoded_bits);
    title('Decoded bits');

    %Compute the actual BER and the theoretical BER for this SNR
    diff=sum(abs(decoded_bits-input));
    BER(d,k)=diff/length(input)
    
    theoretical_awgn_only_BER(k)=qfunc(sqrt(SNR_lin));
    
    %Compute theoretical BER for BPSK modulation with Flat fading channel + AWGN + equalizer
    %if (M==2)
    theoretical_flat_fading_BER(k)=0;
    G=SNR_lin;
    theoretical_flat_fading_BER(k)=(1/2)*(1-sqrt(G/(1+G))); 
    %endif
    
    num_of_transmitter_packets=0;
    
    
    temp_goodput=0;
    for i=1:T*Rb:length(decoded_bits)
      if sum(abs(decoded_bits(i:i+T*Rb-1)-input(i:i+T*Rb-1)))==0
        temp_goodput++;
      endif
      num_of_transmitter_packets++;
      
    endfor
    
    temp_goodput=temp_goodput/num_of_transmitter_packets;
    if MIMO_on==1
      temp_goodput*=MIMO_antennas;
    endif
    
    goodput=[goodput; temp_goodput]
  
  endfor
endfor



figure(8);
semilogy(SNR_db,theoretical_awgn_only_BER,'r','linewidth',2)
hold on;
semilogy(SNR_db,theoretical_flat_fading_BER,'g','linewidth',2)
semilogy(SNR_db,BER(1,:),'b','linewidth',2)
hold off;
legend('AWGN Theoretical BER', 'Flat Fading theoretical BER','Simulation BER');  
   