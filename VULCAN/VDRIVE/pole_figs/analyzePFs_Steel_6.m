%% Import Script for PoleFigure Data

% created from reduce_fit python script @ 2022-07-30 16:40:43.549629
% for MTEX version 5.30

warning('off', 'all')
clear

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
rot2 = rotation('axis',yvector,'angle',0*degree);
rot3 = rotation('axis',xvector,'angle',0*degree);
rot = rot3 * rot2 * rot1;

%% Specify File Names

% path to files
fname_BCC = {...
  'Steel_6_BCC_400.jul',...
  'Steel_6_BCC_321.jul',...
  'Steel_6_BCC_222.jul',...
  'Steel_6_BCC_310.jul',...
  'Steel_6_BCC_220.jul',...
  'Steel_6_BCC_211.jul',...
  'Steel_6_BCC_200.jul',...
  'Steel_6_BCC_110.jul',...
  };

fname_FCC = {...
  'Steel_6_FCC_422.jul',...
  'Steel_6_FCC_420.jul',...
  'Steel_6_FCC_331.jul',...
  'Steel_6_FCC_400.jul',...
  'Steel_6_FCC_222.jul',...
  'Steel_6_FCC_311.jul',...
  'Steel_6_FCC_220.jul',...
  'Steel_6_FCC_200.jul',...
  'Steel_6_FCC_111.jul',...
  };


%% Specify Miller Indice

h_BCC = {...
  Miller(4,0,0,CS_BCC),...
  Miller(3,2,1,CS_BCC),...
  Miller(2,2,2,CS_BCC),...
  Miller(3,1,0,CS_BCC),...
  Miller(2,2,0,CS_BCC),...
  Miller(2,1,1,CS_BCC),...
  Miller(2,0,0,CS_BCC),...
  Miller(1,1,0,CS_BCC),...
  };

h_FCC = {...
  Miller(4,2,2,CS_FCC),...
  Miller(4,2,0,CS_FCC),...
  Miller(3,3,1,CS_FCC),...
  Miller(4,0,0,CS_FCC),...
  Miller(2,2,2,CS_FCC),...
  Miller(3,1,1,CS_FCC),...
  Miller(2,2,0,CS_FCC),...
  Miller(2,0,0,CS_FCC),...
  Miller(1,1,1,CS_FCC),...
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
    allI1{fi} = c.data(:,3);
    allIesd1{fi} = c.data(:,4);
    allIesd1{fi}(allIesd1{fi} > 1E5) = 0;
    allR1{fi} = vector3d.byPolar(deg2rad(c.data(:,1)), deg2rad(c.data(:,2)), 'antipodal');
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
    allI2{fi} = c.data(:,3);
    allIesd2{fi} = c.data(:,4);
    allIesd2{fi}(allIesd2{fi} > 1E5) = 0;
    allR2{fi} = vector3d.byPolar(deg2rad(c.data(:,1)), deg2rad(c.data(:,2)), 'antipodal');
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

    % *** Misorientation eval

    n_gamv     = zeros(n_boot,1);
    n_alpv     = zeros(n_boot,1);
    n_betv     = zeros(n_boot,1);
    n_E_211       = zeros(n_boot,1);
    n_E_311       = zeros(n_boot,1);
    n_tiBCC      = zeros(n_boot,1);
    n_tiFCC      = zeros(n_boot,1);
    n_h1       = length(fname_BCC);
    n_h2       = length(fname_FCC);
    I1 = cell(n_boot,1);
    I2 = cell(n_boot,1);
    cen_rots = equispacedSO3Grid(CS_BCC,SS_BCC,'maxAngle',5*degree,'center',cen_rot,'points',n_boot);

    disp('done with setup.')
    disp('Starting misorientation sampling loop...')
    try
        b = ProgressBar(n_boot, ...
        'IsParallel', true, ...
        'Title', 'Steel_6_texred misori' ...
        );
        updateParallel()
        b.setup([],[],[])
    catch
    end
    
    parfor n=1:n_boot
    
        npf_BCC = rotate(rpf_BCC, cen_rots(n));
        solver   = MLSSolver(npf_BCC);
        nodf_BCC = solver.calcODF('silent');
        n_odfer_BCC(n) = solver.e;
    
        vol_alp = 100 * volume(nodf_BCC,alpha,10*degree);
        vol_gam = 100 * volume(nodf_BCC,gamma,10*degree);

        n_alpv(n) = vol_alp;
        n_gamv(n) = vol_gam;
    
        grains      = nodf_BCC.discreteSample(5000);
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

        npf_FCC = rotate(rpf_FCC, cen_rots(n));

        solver   = MLSSolver(npf_FCC);
        nodf_FCC = solver.calcODF('silent');
        n_odfer_FCC(n) = solver.e;
        
        vol_bet = 100 * volume(nodf_FCC,beta,10*degree);
%         fprintf(' FCC volume fractions | beta fibre: %f6\n',vol_bet)
    
        n_betv(n) = vol_bet;
    
        grains      = nodf_FCC.discreteSample(5000);
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

    out = [n_tiFCC, n_tiBCC, n_odfer_FCC, n_odfer_BCC, n_E_311, n_E_211, n_alpv, n_gamv, n_betv];
    T = array2table(out);
    T.Properties.VariableNames = {'texind_FCC','texind_BCC', 'odferr_FCC', 'odferr_BCC','E_311(FCC)','E_211(BCC)','volfrac_alp','volfrac_gam','volfrac_bet'};
    fname = sprintf('bootstrap_res_misori_n%i_trun.csv',n_boot);
    writetable(T,[pwd '/' fname])

end

%% Import the Data

% create a Pole Figure variable containing the data

% _1 is phase: BCC

pf_1 = PoleFigure.load(fname_1,h_1,CS_BCC,SS_1,'interface','juelich');
odf_1 = calcODF(pf_1);
odf_1 = rotate(odf_1,rot);
save('odf_BCC_Steel_6_G2','odf_1');

% write out rotated pole figures
rotPF_1 = rotate(pf_1,rot);
export(rotPF_1,'BCC','degree')

[odf_1,cen_rot_1,~,~] = centerSpecimen(odf_1,xvector,'Fourier');

figure
plotSection(odf_1,'phi2',45*degree)
xticks([0 45 90 135 180 225 270 315 360])
setColorRange([0 14],'current')
mtexColorbar('title','BCC')
saveFigure('BCC_phi2_45.jpg')
close

% _2 is phase: FCC

pf_2 = PoleFigure.load(fname_2,h_2,CS_FCC,SS_2,'interface','juelich');
odf_2 = calcODF(pf_2);
odf_2 = rotate(odf_2,rot);
save('odf_FCC_Steel_6_G2','odf_2');

% write out rotated pole figures
rotPF_2 = rotate(pf_2,rot);
export(rotPF_2,'FCC','degree')

[odf_2,cen_rot_2,~,~] = centerSpecimen(odf_2,xvector,'Fourier');

figure
plotSection(odf_2,'phi2',45*degree)
xticks([0 45 90 135 180 225 270 315 360])
setColorRange([0 5],'current')
mtexColorbar('title','FCC')
saveFigure('FCC_phi2_45.jpg')
close

%% Plot recalculated pfs

% _1 is phase: BCC

nh_1 = length(h_1);
repf_1 = calcPoleFigure(odf_1,h_1);
figure
plot(repf_1{nh_1-1:nh_1},'contourf',0:0.5:6);
mtexColorbar('title','BCC')
saveFigure('BCC_pfs.jpg')
close

figure
plot(rotPF_1{nh_1})
saveFigure('BCC_rawpfs.jpg')
close

% _2 is phase: FCC

nh_2 = length(h_2);
repf_2 = calcPoleFigure(odf_2,h_2);
figure
plot(repf_2{nh_2-1:nh_2},'contourf',0:0.5:6);
mtexColorbar('title','FCC')
saveFigure('FCC_pfs.jpg')
close

figure
plot(rotPF_2{nh_2})
saveFigure('FCC_rawpfs.jpg')
close

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

vol_alp = 100 * volume(odf_1,alpha,10*degree);
vol_gam = 100 * volume(odf_1,gamma,10*degree);
fprintf(' BCC volume fractions | alpha fibre: %f6, gamma fibre: %f6\n',vol_alp,vol_gam)

vol_bet = 100 * volume(odf_2,beta,10*degree);
fprintf(' FCC volume fractions | beta fibre: %f6\n',vol_bet)

%% texture strength (index) calculation

ti_1 = textureindex(odf_1);
ti_2 = textureindex(odf_2);

fprintf(' Texture index | BCC: %f6 | FCC: %f6 \n',ti_1,ti_2);

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

fodf_1 = FourierODF(odf_1,4);
[TVoigt_1, TReuss_1, THill_1] = calcTensor(fodf_1,T_ijkl);
E_1  = THill_1.YoungsModulus;

figure
plot(E_1,'contourf','complete','upper','colorrange',[115 145]);
mtexColorbar
saveFigure('BCC_betabrassE.jpg')
close

figure
plotPDF(fodf_1,h_1{end},'contourf','colorrange',[0 2])
mtexColorbar
saveFigure('BCC_fodf4_pfs.jpg')
close

fodf_2 = FourierODF(odf_2,4);
[TVoigt_2, TReuss_2, THill_2] = calcTensor(fodf_2,T_ijkl);
E_2  = THill_2.YoungsModulus;

figure
plot(E_2,'contourf','complete','upper','colorrange',[115 145]);
mtexColorbar
saveFigure('FCC_betabrassE.jpg')
close

figure
plotPDF(fodf_2,h_2{end},'contourf','colorrange',[0 2])
mtexColorbar
saveFigure('FCC_fodf4_pfs.jpg')
close

