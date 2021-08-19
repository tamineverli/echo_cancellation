 clc; clear all;	close all;

clear all
[arrayMusica,Fs] = audioread('female_src_1.wav');

sound(arrayMusica,Fs);
count =0;
sizeMusica = size(arrayMusica);

%objCaptado= wavrecord(sizeMusica, Fs);
objCaptado= audiorecorder(Fs,8,1);
disp('Start speaking.')
recordblocking(objCaptado, 10);
disp('End of Recording.');
arrayCaptado = getaudiodata(objCaptado);
%sound(arrayCaptado,Fs);
sizeCaptado=size(arrayCaptado);

filename = 'femaleCaptadoSalaVermelha.wav';
audiowrite(filename,arrayCaptado,Fs);