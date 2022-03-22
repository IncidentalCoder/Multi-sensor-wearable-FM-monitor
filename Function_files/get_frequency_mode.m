function [spectrum,f_vector,f_mode] = get_frequency_mode(data,Fs,trim_length)

% GET_FREQUNCY_MODE Summary of this function goes here
%   This function retuns the main frequency mode of the input data
%   Input variables: data- a vector representing the data that will be analysed
%                    Fs- a scalar value representing the sampling frequency
%                    trim_length- length of triming in Hz. This is used to
%                                 get the frequency mode above a certain frequncy.
%   Output variable: spectrum- A vector of Fourier coefficient values
%                    f_vector- Corresponding frequency vector
%                    f_mode-   a scalar value representing the main frequency mode
%


% ----------- Generating the frequency vector of the spectrum -------------
L = length(data); % data is a column vector
freq_trim_indx = ceil(trim_length*L/Fs+1); % Index upto which trimming needs to be done.
%                                              L/Fs gives number of points/Hz in the frequency spectrum.
%                                              When multiplied by trim_length, it gives the index
%                                              upto which trimming is necessary 
freq_vector = Fs*(0:floor(L/2))/L; % There are L/2 points in the single-sided spectrum.
%                                    Each point will be Fs/L apart.
freq_vector_trimed = freq_vector(freq_trim_indx:end); % frequency vector after trimming

% ----------------------------- Fourier transform -------------------------
FT_2_sided = abs(fft(data))/L; % Two-sided Fourier spectrum. Normalizing by L is generally performed during fft so that it is not neede for inverse fft
FT_data_1_sided = FT_2_sided(1:floor(L/2)+1); % Floor is used for the cases where L is an odd number
FT_data_1_sided(2:end-1) = 2 * FT_data_1_sided(2:end-1); % multiplication by 2 is used to maintain the conservation of energy
FT_data_1_sided_trimed = FT_data_1_sided(freq_trim_indx:end);

% -------------------- Main frequency mode --------------------------------
[~,I] = max(FT_data_1_sided_trimed); % I is the corresponding index
freq_mode_data = freq_vector_trimed(I); % cell matrix with main frequency mode

% ------------------- Finalizing the output variables ---------------------
f_vector = freq_vector_trimed;
spectrum = FT_data_1_sided_trimed;
f_mode = freq_mode_data;
%

end

