%% Initiate Program
clc, clearvars -except inp, close all, f = filesep;

% add path to LIF analysis code
[gitdir] = extractBefore(matlab.desktop.editor.getActiveFilename,'spectrum-model');
addpath(gitdir)
folders = {'sub-programs','DIH','spectrum-model','422-profile','quadrupole-field-model','export-fig','andor-camera-baseline','plotting'};
for i = 1:length(folders)
    addpath(genpath([gitdir f folders{i}]))
end
c = defineConstants();

% data spreadsheet interface
ind = 1;
if ind ~= 1
    %% CLEAR DIRECTORIES
    inp = {};
    openvar('inp')
end
fields = {'date','phase','pol','n','sigx','sigy','Te','mag','t','tau','xc','yc','xL','xR','yL','yR','tE','pwrToInt','shift422','spinpolL','spinpolR','gamL','comment','dir'};
s = cell2struct(inp,fields,2);

model = 're'; % 'fgr' or 're' spectrum model
weightFit = true;
andir = 'an-allpol-re-Ixy-Ic1.7-10.14.21'; % (string) analysis subfolder name

% load spectra
for i = 1:length(s)
    s(i).andir = [s(i).dir f andir];
    load([s(i).andir f 'os.mat'])
    s(i).os = os;
end

maindir = '/mnt/plasma-backup/projects/2021/10.14.21-Polish-LIF-Code';
mkdir(maindir) % make sure main directory exists
figdir = [maindir f 'figs-fits']; % full path for saving spectrum with fits
mkdir(figdir) % make sure figdir exists as a path
savename = 'os-allpol-re.mat';   % file name for fit results structure
linind = find(strcmp({s.pol},'lin'));   % indices of s that correspond to linear imaging polarization
delays = unique([s.t]); % unique delay times in spin polarization data set (units: us)

iter1 = 1:length(delays);
iter2 = 2500:length(s(1).os.local);

load([gitdir f 'spin-pol' f 'p6.17.5-Te40-t100us.mat'],'spinpol')


%% Fit Spectrum with Constrained Model for All LIF-Laser Polarizations

ind = 0;
for i = iter1 % time index
    indt = find([s.t] == delays(i));
    for j = iter2 % local index
            %% Fit Local LIF Spectrum with RE Model
            ind = ind +1;

            dets = {s(indt(1)).os.local(j).dets; s(indt(2)).os.local(j).dets; s(indt(3)).os.local(j).dets};
            spec = {s(indt(1)).os.local(j).spec; s(indt(2)).os.local(j).spec; s(indt(3)).os.local(j).spec};
            specse = {s(indt(1)).os.local(j).specSE; s(indt(2)).os.local(j).specSE; s(indt(3)).os.local(j).specSE};
            x0 = s(linind(i)).os.local(j).x;
            y0 = s(linind(i)).os.local(j).y;
            pol = {s(indt).pol}';
            p0 = spinpol(x0,y0);
            P = p0 - (1-p0);
            if isnan(p0), p0 = 0.8; end
            I = [s(indt(1)).os.batch.I(x0,y0); s(indt(2)).os.batch.I(x0,y0); s(indt(3)).os.batch.I(x0,y0)];
            n = abs(s(linind(i)).os.local(j).n);
            Ge = getGamma(n*1e14,s(i).Te);
            Ti0 = getPohlIonTempForGivenGeAndN(n*1e14,Ge,s(linind(i)).t*1e-6);
            tE = [s(indt(1)).os.tE; s(indt(2)).os.tE; s(indt(3)).os.tE]*1e-9*c.gam422;
            mag = s(linind(i)).os.isPlasmaMagnetized;
            gamL = s(linind(i)).gamL;

            [fit] = fitModelToSpectrumAllPol(dets,spec,specse,x0,y0,pol,Ti0,P,tE,I,model,mag,gamL,weightFit);

            polind = 0;
            names = fieldnames(fit);
            for m = indt
                polind = polind + 1;
                for l = 1:length(names)
                    s(m).os.local(j).fit.(names{l}) = fit(polind).(names{l});
                end
            end

            %% Subplot Fit Results for Quality Control
            fig = figure('Visible','off');
            fig.Position = [179         239        1272         450];

            axind = 0;
            for m = indt
                axind = axind + 1;
                ax{axind} = subplot(1,3,axind);
                dets = s(m).os.local(j).dets;
                spec = s(m).os.local(j).spec;
                specse = s(m).os.local(j).specSE;
                specFit = s(m).os.local(j).fit.spec;
                dEz = s(m).os.local(j).fit.linec;
                a = [min(specFit) max(specFit)/5];
                pol = s(m).pol;

                errorbar(dets,spec,specse,'.-','LineWidth',1.2,'MarkerSize',15)
                hold on
                plot(dets,specFit,'.-','LineWidth',1.2,'MarkerSize',15)
                hold on
                [~,ind5] = sort(dEz);
                for eind = ind5
                    plot([dEz(eind) dEz(eind)],a,'LineWidth',1)
                end

                xlabel('\Delta (MHz)')
                ylabel('LIF Spectrum (a.u.)')
                if isfield(s(m).os.local(j).fit,'eps')
                    eps = num2str(round(sqrt(s(m).os.local(j).fit.eps.*conj(s(m).os.local(j).fit.eps)),2,'significant'));
                    title(['\epsilon = [' eps ']'])
                end

                ax{axind}.Title.Units = 'normalized';
                ax{axind}.Title.Position(2) = ax{axind}.Title.Position(2).*1.05; 
                ax{axind}.PlotBoxAspectRatio = [1 1 1];
                ax{axind}.Position(2) = ax{axind}.Position(2) - .075;

                v(m) = s(m).os.local(j).fit.v;

                txt{m} = text(0,0,s(m).pol);
                txt{m}.Units = 'normalized';
                txt{m}.Position = [0.0238    0.9451         0];
                txt{m}.FontWeight = 'bold';
            end


            an = annotation('textbox');
            an.Position = [0.3809    0.8870    0.2428    0.0638];

            t = s(m).t*1e3;
            p = num2str(s(m).os.local(j).fit.P);
            Ti = num2str(s(m).os.local(j).fit.Ti,3);

            str1 = ['time = ' num2str(t*1e-3) '\mus, x = ' num2str(x0) ' mm, y = ' num2str(y0) ' mm'];
            str2 = ['p_0 = ' p ', T_i = ' Ti ' K'];
            an.String = {[str1 ', ' str2]};
            an.FitBoxToText = 'on';
            an.LineStyle = 'none';
            an.FontSize = 13;
            an.HorizontalAlignment = 'center';

            an2 = annotation('textbox');
            an2.Position = [0.0285    0.0362    0.0589    0.0978];
            str1 = ['i = ' num2str(i)];
            str2 = ['j = ' num2str(j)];
            an2.String = {str1,str2};
            lgd = legend({'data','fit','\sigma-','pi','pi','\sigma+'});
            lgd.Position = [0.9157    0.5398    0.0665    0.2281];

            an3 = annotation('textbox');
            an3.Position = [0.0168    0.8133    0.0939    0.1400];
            str1 = ['v1 = ' num2str(round(v(indt(1)),2,'significant')) ' m/s'];
            str2 = ['v2 = ' num2str(round(v(indt(2)),2,'significant')) ' m/s'];
            str3 = ['v3 = ' num2str(round(v(indt(3)),2,'significant')) ' m/s'];
            an3.String = {str1,str2,str3};


            lgd = legend({'data','fit','\sigma-','pi','pi','\sigma+'});
            lgd.Position = [0.9157    0.5398    0.0665    0.2281];

            saveas(fig,[figdir f num2str(s(m).t) 'us_' num2str(j) '.png'])
%             save([maindir f savename],'s','-mat')
            close(fig)
            %% Output analysis progress

            disp(['Complete: ' num2str(round(ind/length(iter1)/length(iter2)*100,3,'significant'))  '%'])
    end
    
    % save fit results after each time point
    for j = indt
        os = s(j).os;
        save([s(j).andir f savename],'os','-mat')
    end
end

%% Process Fit Results
for i = 1:length(s)
    for j = 1:length(s(i).os.local)
        x0 = s(i).os.local(j).x;
        y0 = s(i).os.local(j).y;
        
        % Get density from fit by getting fit-amplitude-to-density conversion factor using the
        % central region
        dets = s(i).os.local(j).dets;
        spec = s(i).os.local(j).spec - s(i).os.local(j).fit.offset;
        specfit = s(i).os.local(j).fit.spec - s(i).os.local(j).fit.offset;
        intspec = trapz(dets,spec);
        intspecfit = trapz(dets,specfit);
        
        s(i).os.local(j).nLIF1 = intspec*s(i).os.imgs.intSigToDens(x0,y0);
        s(i).os.local(j).nLIF2 = intspecfit*s(i).os.imgs.intSigToDens(x0,y0);
                
        s(i).os.local(j).nFit = s(i).os.convFacAmpToDens*s(i).os.local(j).fit.amp;
        s(i).os.local(j).nFitErr = s(i).os.convFacAmpToDens*s(i).os.local(j).fit.ampErr;
        
        s(i).os.local(j).intspec = intspec;
        s(i).os.local(j).intspecfit = intspecfit;
        
        s(i).os.local(j).vExp = s(i).os.local(j).fit.v;
        s(i).os.local(j).Ti = s(i).os.local(j).fit.Ti;
        s(i).os.local(j).P = s(i).os.local(j).fit.P;
    end
    
    fieldstr = {'nLIF1','nLIF2','nFit','intspec','intspecfit','vExp','Ti','P'};
    map.x = sort(unique([s(i).os.local.x]));
    map.y = sort(unique([s(i).os.local.y]));
    [map.X,map.Y] = meshgrid(map.x,map.y);

    for j = 1:length(fieldstr)
        map.(fieldstr{j}) = list2mat([s(i).os.local.x],[s(i).os.local.y],[s(i).os.local.(fieldstr{j})],map.X,map.Y);
    end
    
    s(i).os.map = map;
    
    os = s(i).os;
    save([s(i).andir f savename],'os','-mat')
end


