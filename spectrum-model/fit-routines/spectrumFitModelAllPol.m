function [spec] = spectrumFitModelAllPol(dets,ind,amp,o1,o2,o3,v1,v2,v3,Ti,P,I,pol,B,theta,phi,tE,radpat,mod,gamL,gamD,c)
o = [o1 o2 o3]; % offset for each polarization
v = [v1 v2 v3]; % velocity guess for each polarization

spec = zeros(size(dets));
for i = 1:length(pol)
    s = spectrumModel(dets(i == ind),v(i),Ti,P,I(i),pol{i},B,theta,phi,tE(i),radpat,mod,gamL,gamD,c);
    spec(i == ind) = [s.S].*amp+o(i);
end

end