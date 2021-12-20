function [spec] = spectrumFitModel(dets,amp,offset,v,Ti,P,I,pol,B,theta,phi,tE,radpat,mod,gamL,gamD,c)
% this function sets up the spectrum model to have an amplitude and offset for fitting purposes
% see 'spectrumModel.m' for model specifics
s = spectrumModel(dets,v,Ti,P,I,pol,B,theta,phi,tE,radpat,mod,gamL,gamD,c);
spec = s.S.*amp + offset;
end