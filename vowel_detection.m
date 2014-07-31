%Recognition of Speech 
%Scify Prject
%Think Freedom Project

clc;clear;
           load mtlb;

           filename = 'sentence.wav';
   
           Fs = 7418;   % sampling freq (hz)
           duration =3; % duration (sec)
%          z = wavrecord(duration*Fs,Fs);
%          wavwrite(z,Fs,filename);
%          x= wavread('a.wav');
%          x=mtlb;
%This loop is used to take the average frequency of user's frequency.

%So that to know the number of loops
nof_vowels=input('Please enter number of vowels that you want to recognize: ');
names = cell(nof_vowels,1);
for k=1:nof_vowels
name_of_vowel=input('Please enter the name of the vowel: ','s')
%we store in a array of strings the name of the vowels


names{k} = name_of_vowel;

x= wavread(name_of_vowel);



%if i dont know the duration of a wav file i can  compute it         
%            TotalTime = length(x) ./ FS 
%Use the spectrogram to identify a voiced segment for analysis.  
segmentlen = 100;
noverlap = 90;
NFFT = 128;
[y,f,t,p] = spectrogram(x,segmentlen,noverlap,NFFT,Fs);
surf(t,f,10*log10(abs(p)),'EdgeColor','none');
axis xy; axis tight; colormap(jet); view(0,90);
xlabel('Time');
ylabel('Frequency (Hz)');
% we initilize the variables of formants 1 and 2
formants1=0;
formants2=0;
%it counts the number of iteration
counter=0;
%The variable first_ is showing where we start to take th sample
first_s=0.1; 
%We take a sample every 0.25 to analyze the whole .wav file
for first_s=first_s:0.25:duration
    %here we compute where the sample ends 
    
        next_s=first_s+0.15;
%     wecomute the number of iteration
        counter= counter+1;
%     if next_s>duration
%        next_s=next_s-0.05;
%        display('TRIM')
%     end
    
    
%Extract the segment from 0.1 to 0.25 seconds for analysis.
dt = 1/Fs;
I0 = round(first_s/dt);
Iend = round(next_s/dt);
xx = x(I0:Iend);


%Window the speech segment using a Hamming window.


x1 = xx.*hamming(length(xx));


%Apply a pre-emphasis filter. The pre-emphasis filter is a highpass all-pole (AR(1)) filter.
preemph = [1 0.63];
x1 = filter(1,preemph,x1);

%Obtain the linear prediction coefficients.
A = lpc(x1,8);
rts = roots(A);

% Retain only the roots with one sign for the imaginary part and determine the angles corresponding to the roots.
rts = rts(imag(rts)>=0);
angz = atan2(imag(rts),real(rts));

%Convert the angular frequencies in radians/sample represented by the angles to hertz and calculate the bandwidths of the formants.
[frqs,indices] = sort(angz.*(Fs/(2*pi)));
bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)));

%Use the criteria that formant frequencies should be greater than 90 Hz with bandwidths less than 400 Hz to determine the formants.
nn = 1;
for kk = 1:length(frqs)
    if (frqs(kk) > 90 && bw(kk) <400)
        formants(nn) = frqs(kk);
        nn = nn+1;
        
    end
end

formants1=formants1+formants(1);
formants2=formants2+formants(2);



end

% disp('Formant 1 is ')
%  formants1
% disp('Formant 2 is ') 
%  formants2
 
 formants1=formants1/counter;
 formants2=formants2/counter;

%in the struct we store tthe formant 1 and formant2 for every vowel 
  struct.names{k}.final_formant1=formants1;
  struct.names{k}.final_formant2=formants2;

end

duration_sentence=input('Please enter the duration of your sentence: ');

disp('Please press enter to record your sentence:');
pause;
%We record the new sentence that we want to analyze
% z=wavrecord(duration_sentence*Fs,Fs);
% wavwrite(z,Fs,filename);
%We analyze the sentence
z=wavread(filename);
plot(z)
segmentlen = 100;
noverlap = 90;
NFFT = 128;
[y,f,t,p] = spectrogram(z,segmentlen,noverlap,NFFT,Fs);
surf(t,f,10*log10(abs(p)),'EdgeColor','none');
axis xy; axis tight; colormap(jet); view(0,90);
xlabel('Time');
ylabel('Frequency (Hz)');
%The variable first_ is showing where we start to take th sample
first_s=0.1; 
%We take a sample every 0.25 to analyze the whole .wav file
for first_s=first_s:0.25:duration_sentence
    %here we compute where the sample ends 
    
    next_s=first_s+0.15;
%     if next_s>duration
%        next_s=next_s-0.05;
%        display('TRIM')
%     end
    
    
%Extract the segment from 0.1 to 0.25 seconds for analysis.
dt = 1/Fs;
I0 = round(first_s/dt);
Iend = round(next_s/dt);
xx = z(I0:Iend);


%Window the speech segment using a Hamming window.


x1 = xx.*hamming(length(xx));


%Apply a pre-emphasis filter. The pre-emphasis filter is a highpass all-pole (AR(1)) filter.
preemph = [1 0.63];
x1 = filter(1,preemph,x1);

%Obtain the linear prediction coefficients.
A = lpc(x1,8);
rts = roots(A);

% Retain only the roots with one sign for the imaginary part and determine the angles corresponding to the roots.
rts = rts(imag(rts)>=0);
angz = atan2(imag(rts),real(rts));

%Convert the angular frequencies in radians/sample represented by the angles to hertz and calculate the bandwidths of the formants.
[frqs,indices] = sort(angz.*(Fs/(2*pi)));
bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)));

%Use the criteria that formant frequencies should be greater than 90 Hz with bandwidths less than 400 Hz to determine the formants.
nn = 1;
for kk = 1:length(frqs)
    if (frqs(kk) > 90 && bw(kk) <400)
        formants(nn) = frqs(kk);
        nn = nn+1;
        
    end
end

 %We check for every vowel that we have record ,if the sentence respond to any vowel   
   
    for jj=1:k
%         disp('Formant1 is: ')
     formants(1);
%      disp('Formant2 is: ')
     formants(2);
    if (((struct.names{jj}.final_formant1 - 40)<formants(1) && (struct.names{jj}.final_formant1 + 40)>formants(1)) && ((struct.names{jj}.final_formant2 - 40)<formants(2) && (struct.names{jj}.final_formant2 + 40)>formants(2)))
        
           S=sprintf('The vowel that was detected is  %s',names{jj});
           disp(S)
     
   
    end
    
    end
end



clear y F;