%% Import Script for PoleFigure Data

% created from reduce_fit python script @ 2022-07-13 17:21:00.264875
% for MTEX version 5.30

warning('off', 'all')
clear
clc

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS_BCC = crystalSymmetry('m-3m', [1 1 1]); % lattice-type:I
SS_BCC = specimenSymmetry('1');
CS_FCC = crystalSymmetry('m-3m', [1 1 1]); % lattice-type:F
SS_FCC = specimenSymmetry('1');

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

rot1 = rotation('axis',zvector,'angle',0*degree);
rot2 = rotation('axis',yvector,'angle',90*degree);
rot = rot2 * rot1;

%% Specify File Names

% path to files
fname_BCC = {...
  'Steel_6_texred_uq_BCC_220.pf',...
  'Steel_6_texred_uq_BCC_200.pf',...
  };

fname_FCC = {...
  'Steel_6_texred_uq_FCC_222.pf',...
  'Steel_6_texred_uq_FCC_311.pf',...
  'Steel_6_texred_uq_FCC_220.pf',...
  };


%% Specify Miller Indice

h_BCC = {...
  Miller(2,2,0,CS_BCC),...
  Miller(2,0,0,CS_BCC),...
  };

h_FCC = {...
  Miller(2,2,2,CS_FCC),...
  Miller(3,1,1,CS_FCC),...
  Miller(2,2,0,CS_FCC),...
  };

%% Volume fractions of texture components

BFOri = [[35.3 45.0 0.0 ];...
           [33.6 47.7 5.0 ];...
            [32.1 51.0 10.0];...
            [31.1 54.7 15.0];...
            [31.3 59.1 20.0];...
            [35.2 64.2 25.0];...
            [46.1 69.9 30.0];...
            [49.5 76.2 35.0];...
            [51.8 83.0 40.0];...
            [54.7 90.0 45.0];...
            [90.0 35.3 45.0];...
            [80.2 35.4 50.0];...
            [73.0 35.7 55.0];...
            [66.9 36.2 60.0];...
            [61.2 37.0 65.0];...
            [55.9 38.0 70.0];...
            [50.7 39.2 75.0];...
            [45.6 40.8 80.0];...
            [40.5 42.7 85.0];...
            [35.3 45.0 90.0]]*degree;

beta = orientation('Euler', BFOri(:,1), BFOri(:,2), BFOri(:,3),crystalSymmetry('m-3m'),specimenSymmetry('1'),'bunge');
alpha = fibre.alpha(crystalSymmetry('m-3m'),specimenSymmetry('1'),'full');
gamma = fibre.gamma(crystalSymmetry('m-3m'),specimenSymmetry('1'),'full');

%% Import the Data

for fi=1:length(fname_BCC)
    c = importdata(fname_BCC{fi}, ' ', 3);
    allI1{fi} = c.data(:,4);
    allIesd1{fi} = c.data(:,5);
    allIesd1{fi}(allIesd1{fi} > 1E5) = 0;
    allR1{fi} = vector3d(c.data(:,1), c.data(:,2), c.data(:,3), 'antipodal');
    allH1{fi} = h_BCC{fi};
end

pf_BCC = PoleFigure(allH1, allR1, allI1);
rpf_BCC = rotate(pf_BCC, rot);
odf_BCC = calcODF(pf_BCC);
rodf_BCC = rotate(odf_BCC,rot);
repf_BCC = calcPoleFigure(rodf_BCC, h_BCC);
[~,cen_rot,~,~] = centerSpecimen(rodf_BCC,xvector,'Fourier','silent');


for fi=1:length(fname_FCC)
    c = importdata(fname_FCC{fi}, ' ', 3);
    allI2{fi} = c.data(:,4);
    allIesd2{fi} = c.data(:,5);
    allIesd2{fi}(allIesd2{fi} > 1E5) = 0;
    allR2{fi} = vector3d(c.data(:,1), c.data(:,2), c.data(:,3), 'antipodal');
    allH2{fi} = h_FCC{fi};
end

pf_FCC = PoleFigure(allH2, allR2, allI2);
rpf_FCC = rotate(pf_FCC, rot);
odf_FCC = calcODF(pf_FCC);
rodf_FCC = rotate(odf_FCC,rot);
repf_FCC = calcPoleFigure(rodf_FCC, h_FCC);

%% betaBrass elasticity homogenization 
C11 = 204.6;
C12 = 137.7;
C13 = C12;
C23 = C12;
C22 = C11;
C33 = C22;
C44 = 126.2;
C55 = C44;
C66 = C44;
T_cs = crystalSymmetry('m-3m', [1 1 1], 'mineral',...
    'Copper', 'color', 'light blue');

c_ijkl = stiffnessTensor(...
    [[  C11     C12    C13    0.0     0.0    0.0];...
    [   C12     C22    C23    0.0     0.0    0.0];...
    [   C13     C23    C33    0.0     0.0    0.0];...
    [   0.0     0.0    0.0    C44     0.0    0.0];...
    [   0.0     0.0    0.0    0.0     C55    0.0];...
    [   0.0     0.0    0.0    0.0     0.0    C66]],T_cs);

s_ijkl = inv(c_ijkl);

n_boots = [5000];

%setup the pool
try
    parpool;
catch
end

for bi=1:length(n_boots)

    n_boot = n_boots(bi);

    ncen_rot_BCC = cell(n_boot,1);
    n_gamv     = zeros(n_boot,1);
    n_alpv     = zeros(n_boot,1);
    n_betv     = zeros(n_boot,1);
    n_odfer_FCC   = zeros(n_boot,1);
    n_odfer_BCC   = zeros(n_boot,1);
    n_E_211       = zeros(n_boot,1);
    n_E_311       = zeros(n_boot,1);
    n_tiBCC      = zeros(n_boot,1);
    n_tiFCC      = zeros(n_boot,1);
    n_h1       = length(fname_BCC);
    n_h2       = length(fname_FCC);
    I1 = cell(n_boot,1);
    I2 = cell(n_boot,1);

    for n=1:n_boot

        I1{n} = cell(1,n_h1);
        I2{n} = cell(1,n_h2);

        % odf_BCC
        for fi=1:n_h1
            I1{n}{fi} = normrnd(allI1{fi}, allIesd1{fi});
        end

        % odf_FCC
        for fi=1:n_h2
            I2{n}{fi} = normrnd(allI2{fi}, allIesd2{fi});
        end

    end

    disp('done with setup.')
    disp('Starting peak intensity sampling loop...')
    
    try
        b = ProgressBar(n_boot, ...
        'IsParallel', true, ...
        'Title', 'Steel_6_texred' ...
        );
        updateParallel()
    
        b.setup([],[],[])
    catch
    end

    parfor n=1:n_boot
    
        npf_BCC  = PoleFigure(allH1, allR1, I1{n});
%         npf_BCC(npf_BCC.isOutlier) = [];

        solver   = MLSSolver(npf_BCC);
        nodf_BCC = solver.calcODF('silent');
        nodf_BCC = rotate(nodf_BCC,rot);
        [~,ncen_rot_BCC{n},~,~] = centerSpecimen(nodf_BCC,xvector,'Fourier','silent');
        nodf_BCC = rotate(nodf_BCC, ncen_rot_BCC{n});
    
        vol_alp = 100 * volume(nodf_BCC,alpha,10*degree);
        vol_gam = 100 * volume(nodf_BCC,gamma,10*degree);

        n_alpv(n) = vol_alp;
        n_gamv(n) = vol_gam;
    
        grains      = nodf_BCC.discreteSample(5000);
        inv_g       = inv(grains);
        THill_HEM   = calcTensor(grains, c_ijkl, 'Hill');
        diff_fibre  = fibre(Miller(2,1,1,CS_BCC), vector3d.Y);
        ind         = angle(diff_fibre, grains) < 10*degree;
        diff_grains = grains(ind); 
        THill_diff  = calcTensor(diff_grains, s_ijkl);
        THill = THill_diff+inv(THill_HEM);
        
        m = [0, 1, 0];
        f = 0;
        for i=1:3
            for j=1:3
                % Shortcut because we only care about 22 stress
                k = 2;
                l = 2;
                f = f + 0.5*m(i)*m(j)*THill.M(i,j,k,l);
            end
        end
    
        n_E_211(n) = 1/f;
    
        npf_FCC  = PoleFigure(allH2, allR2, I2{n});
%         npf_FCC(npf_FCC.isOutlier) = [];
    
        solver   = MLSSolver(npf_FCC);
        nodf_FCC = solver.calcODF('silent');

        nodf_FCC = rotate(nodf_FCC,rot);
        % Use the centering from odf_BCC
        nodf_FCC = rotate(nodf_FCC,ncen_rot_BCC{n});
    
        vol_bet = 100 * volume(nodf_FCC,beta,10*degree);

    
        n_betv(n) = vol_bet;
    
        grains      = nodf_FCC.discreteSample(5000);
        inv_g       = inv(grains);
        THill_HEM   = calcTensor(grains, c_ijkl, 'Hill');
        diff_fibre  = fibre(Miller(3,1,1,CS_FCC), vector3d.Y);
        ind         = angle(diff_fibre, grains) < 10*degree;
        diff_grains = grains(ind); 
        THill_diff  = calcTensor(diff_grains, s_ijkl);
        THill = THill_diff+inv(THill_HEM);
        
        m = [0, 1, 0];
        f = 0;
        for i=1:3
            for j=1:3
                % Shortcut because we only care about 22 stress
                k = 2;
                l = 2;
                f = f + 0.5*m(i)*m(j)*THill.M(i,j,k,l);
            end
        end
    
        n_E_311(n) = 1/f;
    
        n_tiBCC(n) = textureindex(nodf_BCC);
        n_tiFCC(n) = textureindex(nodf_FCC);
         
        try
            updateParallel();
        catch
        end

    end
    try
        b.release()
    catch
    end

    for i=1:length(ncen_rot_BCC)
        cen_rots(i,:) = ncen_rot_BCC{i}.Euler;
    end
    
    out = [n_tiFCC, n_tiBCC, n_E_311, n_E_211, n_alpv, n_gamv, n_betv, cen_rots];
    T = array2table(out);
    T.Properties.VariableNames = {'texind_FCC','texind_BCC','E_311(FCC)','E_211(BCC)','volfrac_alp','volfrac_gam','volfrac_bet','phi1','Phi','phi2'};
    fname = sprintf('bootstrap_res_n%i_trun.csv',n_boot);
    writetable(T,[pwd '/' fname])

    clear cen_rots

end

%%

% write out the odf
rodf_BCC = rotate(rodf_BCC,cen_rot);
save('odf_BCC_Steel_6_texred_uq','rodf_BCC');

% write out rotated pole figures
rpf_BCC  = rotate(rpf_BCC, cen_rot);
export(rpf_BCC,'BCC','degree')

figure
plotSection(rodf_BCC,'phi2',45*degree)
xticks([0 45 90 135 180 225 270 315 360])
setColorRange([0 14],'current')
mtexColorbar('title','BCC')
saveFigure('BCC_phi2_45.jpg')
close

% _2 is phase: FCC

rodf_FCC = rotate(rodf_FCC,cen_rot);
save('odf_FCC_Steel_6_texred_uq','rodf_FCC');

% write out rotated pole figures
rpf_FCC  = rotate(rpf_FCC, cen_rot);
export(rpf_FCC,'FCC','degree')

figure
plotSection(rodf_FCC,'phi2',45*degree)
xticks([0 45 90 135 180 225 270 315 360])
setColorRange([0 5],'current')
mtexColorbar('title','FCC')
saveFigure('FCC_phi2_45.jpg')
close

ti_BCC = textureindex(rodf_BCC);
ti_FCC = textureindex(rodf_FCC);
fprintf(' Texture index | BCC: %f6 | FCC: %f6 \n',ti_BCC,ti_FCC);

vol_alp = 100 * volume(rodf_BCC,alpha,10*degree);
vol_gam = 100 * volume(rodf_BCC,gamma,10*degree);
fprintf(' BCC volume fractions | alpha fibre: %f6, gamma fibre: %f6\n',vol_alp,vol_gam)

vol_bet = 100 * volume(rodf_FCC,beta,10*degree);
fprintf(' FCC volume fractions | beta fibre: %f6\n',vol_bet)

grains      = rodf_BCC.discreteSample(5000);
inv_g       = inv(grains);
THill_HEM   = calcTensor(grains, c_ijkl, 'Hill');
diff_fibre  = fibre(Miller(2,1,1,CS_BCC), vector3d.Y);
ind         = angle(diff_fibre, grains) < 10*degree;
diff_grains = grains(ind); 
THill_diff  = calcTensor(diff_grains, s_ijkl);
THill = THill_diff+inv(THill_HEM);

m = [0, 1, 0];
f = 0;
for i=1:3
    for j=1:3
        % Shortcut because we only care about 22 stress
        k = 2;
        l = 2;
        f = f + 0.5*m(i)*m(j)*THill.M(i,j,k,l);
    end
end

E_211 = 1/f;

grains      = rodf_FCC.discreteSample(5000);
inv_g       = inv(grains);
THill_HEM   = calcTensor(grains, c_ijkl, 'Hill');
diff_fibre  = fibre(Miller(3,1,1,CS_BCC), vector3d.Y);
ind         = angle(diff_fibre, grains) < 10*degree;
diff_grains = grains(ind); 
THill_diff  = calcTensor(diff_grains, s_ijkl);
THill = THill_diff+inv(THill_HEM);

m = [0, 1, 0];
f = 0;
for i=1:3
    for j=1:3
        % Shortcut because we only care about 22 stress
        k = 2;
        l = 2;
        f = f + 0.5*m(i)*m(j)*THill.M(i,j,k,l);
    end
end

E_311 = 1/f;

out = [ti_FCC, ti_BCC, E_311, E_211, vol_alp, vol_gam, vol_bet];
T = array2table(out);
T.Properties.VariableNames = {'texind_FCC','texind_BCC','E_311(FCC)','E_211(BCC)','volfrac_alp','volfrac_gam','volfrac_bet'};
fname = sprintf('res_trun.csv');
writetable(T,[pwd '/' fname])

%% Plot recalculated pfs

% _1 is phase: BCC

nh_BCC = length(h_BCC);
repf_1 = calcPoleFigure(rodf_BCC,h_BCC);
figure
plot(repf_1{nh_BCC-1:nh_BCC},'contourf',0:0.5:5.5);
mtexColorbar('title','BCC')
saveFigure('BCC_pfs.jpg')
close

figure
plot(rpf_BCC{nh_BCC-1:nh_BCC})
saveFigure('BCC_rawpfs.jpg')
close

% _2 is phase: FCC

nh_FCC = length(h_FCC);
repf_2 = calcPoleFigure(rodf_FCC,h_FCC);
figure
plot(repf_2{nh_FCC-1:nh_FCC},'contourf',0:0.5:4);
mtexColorbar('title','FCC')
saveFigure('FCC_pfs.jpg')
close

figure
plot(rpf_FCC{nh_FCC-1:nh_FCC})
saveFigure('FCC_rawpfs.jpg')
close

%% betaBrass elasticity homogenization 
C11 = 52;
C12 = 34;
C13 = C12;
C23 = C12;
C22 = 52;
C33 = C22;
C44 = 173;
C55 = C44;
C66 = C44;
T_cs = crystalSymmetry('m-3m', [1 1 1], 'mineral',...
    'Copper', 'color', 'light blue');

T_ijkl = stiffnessTensor(...
    [[  C11     C12    C13    0.0     0.0    0.0];...
    [   C12     C22    C23    0.0     0.0    0.0];...
    [   C13     C23    C33    0.0     0.0    0.0];...
    [   0.0     0.0    0.0    C44     0.0    0.0];...
    [   0.0     0.0    0.0    0.0     C55    0.0];...
    [   0.0     0.0    0.0    0.0     0.0    C66]],T_cs);

fodf_1 = FourierODF(rodf_BCC,4);
[TVoigt_1, TReuSS_BCC, THill_1] = calcTensor(fodf_1,T_ijkl);
E_1  = THill_1.YoungsModulus;

figure
plot(E_1,'contourf','complete','upper','colorrange',[115 145]);
mtexColorbar
saveFigure('BCC_betabrassE.jpg')
close

figure
plotPDF(fodf_1,h_BCC{end},'contourf','colorrange',[0 2])
mtexColorbar
saveFigure('BCC_fodf4_pfs.jpg')
close

fodf_2 = FourierODF(rodf_FCC,4);
[TVoigt_2, TReuSS_FCC, THill_2] = calcTensor(fodf_2,T_ijkl);
E_2  = THill_2.YoungsModulus;

figure
plot(E_2,'contourf','complete','upper','colorrange',[115 145]);
mtexColorbar
saveFigure('FCC_betabrassE.jpg')
close

figure
plotPDF(fodf_2,h_FCC{end},'contourf','colorrange',[0 2])
mtexColorbar
saveFigure('FCC_fodf4_pfs.jpg')
close


