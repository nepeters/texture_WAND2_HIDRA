%% Load in MAUD ODFs
% using import tool written by Dan Savage
% https://github.com/djm87/MAUD-TBP

% main directory
parfile_dir = './samples_texred_5deg';
prefix = '_texred_5deg';

CS_bcc = crystalSymmetry('m-3m', [2.9 2.9 2.9]);
CS_fcc = crystalSymmetry('m-3m', [3.6 3.6 3.6]);

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

betaFibre = orientation('Euler', BFOri(:,1), BFOri(:,2), BFOri(:,3),CS_fcc,specimenSymmetry('1'),'bunge');
     
% Specify Miller Indice
h_bcc = { ...
  Miller(1,1,0,CS_bcc),...
  Miller(2,0,0,CS_bcc),...
  Miller(2,1,1,CS_bcc),...
  Miller(2,2,0,CS_bcc),...
  Miller(3,1,0,CS_bcc),...
  Miller(2,2,2,CS_bcc),...
  Miller(3,2,1,CS_bcc),...
  Miller(4,0,0,CS_bcc),...
  };

h_fcc = { ...
  Miller(1,1,1,CS_fcc),...
  Miller(2,0,0,CS_fcc),...
  Miller(2,2,0,CS_fcc),...
  Miller(3,1,1,CS_fcc),...
  Miller(2,2,2,CS_fcc),...
  Miller(4,0,0,CS_fcc),...
  Miller(3,3,1,CS_fcc),...
  Miller(4,2,0,CS_fcc),...
  Miller(4,2,2,CS_fcc),...
  };

h = {h_bcc, h_fcc};

% loop through steel samples
odf = {};
new_rot = {};

% original angles I used.. need to adjust a few
orig_rotAng = {[0 0],[88 5],[90 25],[88 65],[0 73],[90 0],[0 17],[90 0],[90 105],[90 95]};

for i=2:10
    
    sample = ['Steel_' int2str(i)];
    sample_file = [parfile_dir '/' sample '/' 'HB2C_MAUD_' sample '_Step4.par'];
    % make an image folder
    status = mkdir([parfile_dir '/' sample '/Figures']);
    
    disp(['loading ODF(s) from' sample]);
    odf_in = ExtractODFFromPar(sample_file);
    
    % store in cell
    odf{i} = {};
    % loop over phases
    for pi=1:odf_in.num

        maud_odf = odf_in.odf{pi};
        name_odf = ['Steel_' int2str(i) '-' odf_in.name{pi}];
        
        % === uncentered cells ===
        % mask Beta == 0
        mask_bz = (maud_odf(:,2) == 0);
        % mask Beta == max
        mask_bm = (maud_odf(:,2) == max(maud_odf(:,2)));
        
        % mask Alpha == 0
        mask_az = (maud_odf(:,1) == 0);
        % mask Alpha == max
        mask_am = (maud_odf(:,1) == max(maud_odf(:,1)));       

        % mask Gamma == 0
        mask_gz = (maud_odf(:,3) == 0);
        % mask Gamma == max
        mask_gm = (maud_odf(:,3) == max(maud_odf(:,3)));
        
        da_dg = ones(length(maud_odf),1) * deg2rad(odf_in.res(pi))^2;
        
        % alpha edge cases - 0.5*Δalpha * Δgamma
        da_dg(or(mask_az,mask_am)) = 0.5*deg2rad(odf_in.res(pi)) * deg2rad(odf_in.res(pi));
        % gamma edge cases - 0.5*Δalpha * Δgamma
        da_dg(or(mask_gz,mask_gm)) = 0.5*deg2rad(odf_in.res(pi)) * deg2rad(odf_in.res(pi));
        % alpha and gamma edge case - 0.5*Δalpha * 0.5*Δgamma
        da_dg(and(or(mask_gz,mask_gm),or(mask_az,mask_am))) = 0.5*deg2rad(odf_in.res(pi)) * 0.5*deg2rad(odf_in.res(pi)); 

        delta_Phi = ones(length(maud_odf),1) .* ( cos( deg2rad(maud_odf(:,2)) - ( deg2rad(odf_in.res(pi))/2 ) ) - cos( deg2rad(maud_odf(:,2)) + ( deg2rad(odf_in.res(pi))/2 ) ) );
        % Beta = 0
        delta_Phi(mask_bz) = ( cos( deg2rad(maud_odf(mask_bz,2))) - cos( deg2rad(maud_odf(mask_bz,2)) + ( deg2rad(odf_in.res(pi))/2 ) ) );
        % Beta = max
        delta_Phi(mask_bm) = ( cos( deg2rad(maud_odf(mask_bm,2)) - ( deg2rad(odf_in.res(pi))/2 ) ) - cos( deg2rad(maud_odf(mask_bm,2)) ) );       
        cellVolume = da_dg .* delta_Phi;
        
%         % debug plot
%         color = cellVolume*10000;
%         size  = cellVolume*10000;
%         colormap(spring)
%         scatter3(maud_odf(:,1),maud_odf(:,2),maud_odf(:,3),size,color)
%         
%         % === centered cells ===
%         delta_Phi = ones(length(maud_odf),1) .* ( cos( deg2rad(maud_odf(:,2)) - ( deg2rad(odf_in.res(pi))/2 ) ) - cos( deg2rad(maud_odf(:,2)) + ( deg2rad(odf_in.res(pi))/2 ) ) );
%         cellVolume = deg2rad(odf_in.res(pi)) * deg2rad(odf_in.res(pi)) * delta_Phi;

        temp = sum( maud_odf(:,4) .* cellVolume ) / sum( cellVolume );
        orig_wgts = maud_odf(:,4);
        
        if temp ~= 1

            od_dg = maud_odf(:,4) .* cellVolume;
            norm = sum( od_dg ) / sum( cellVolume );
            weights = ( 1 / norm ) .* maud_odf(:,4);
        
        end
        
        maud_odf(:,4) = maud_odf(:,4) .* cellVolume;
            
        odf_out = [parfile_dir '/' sample '/' 'HB2C_MAUD_' sample '_' odf_in.name{pi} '.txt'];
        fout = fopen(odf_out,'w+');
        fprintf(fout,'alpha\tbeta\tgamma\tweight\n');
        fprintf(fout,'%g\t %g\t %g\t %g\n',maud_odf');

        % specify kernel
        psi = deLaValleePoussinKernel('halfwidth',5*degree);

        % load the ODF into the variable odf
        odf{i}{pi} = ODF.load(odf_out,odf_in.cs{pi},odf_in.ss{pi},'density','kernel',psi,'resolution',5*degree,...
          'interface','generic',...
          'ColumnNames', { 'alpha' 'beta' 'gamma' 'weight'}, 'Matthies', 'Degrees');

        odf{i}{pi} = rotate(odf{i}{pi},rotation('axis',zvector,'angle',90*degree));
      
        maud_ori = orientation.byEuler(maud_odf(:,1)*degree,maud_odf(:,2)*degree,maud_odf(:,3)*degree,'ZYZ',odf_in.cs{pi},odf_in.ss{pi});

        values=odf{i}{pi}.eval(maud_ori);
        err=sum(abs(orig_wgts-values));
        [max_err,ind]=max(abs(orig_wgts-values));
        mean_err=mean(abs(orig_wgts-values));
        fprintf(' %s | Import err : %f6, err max: %f6, err mean: %f6\n',name_odf,err,max_err,mean_err)      

    end
    
    if i > 6
        % FCC is majority
        [~,cent_rot,~,~] = centerSpecimen(odf{i}{2},xvector,'Fourier');
    else
        % BCC is majority
        [~,cent_rot,~,~] = centerSpecimen(odf{i}{1},xvector,'Fourier');    
    end

    for pi=1:odf_in.num
    
        maud_odf = odf_in.odf{pi};
        name_odf = ['Steel_' int2str(i) '-' odf_in.name{pi}];
        
        odf{i}{pi} = rotate(odf{i}{pi},cent_rot);

        disp('plotting ODF phi2=45deg section')

        figure
        plotSection(odf{i}{pi},'phi2',45*degree)
        xticks([0 45 90 135 180 225 270 315 360])
        if strcmp('Fe(FCC)',odf_in.name{pi})
            setColorRange([0 5],'current')
        else
            setColorRange([0 14],'current')
        end
        mtexColorbar('title',odf_in.name{pi})
        fname = [parfile_dir '/' sample '/Figures/' name_odf prefix '_phi2_45.jpg'];
        saveFigure(fname)
        close
        
        export(odf{i}{pi},[parfile_dir '/' sample '/Figures/' name_odf prefix '_ODF.MTEX'],...
            'generic','Bunge','resolution',5*degree)

        disp('plotting recalculated pole figures')

        figure

        if pi == 1
            plotPDF(odf{i}{pi},h{pi}(1:3),'contourf',0:0.5:5.5)

            % look at fibre volume fractions
            vol_alp = 100 * volume(odf{i}{pi},fibre.alpha(CS_bcc,specimenSymmetry('1'),'full'),10*degree);
            vol_gam = 100 * volume(odf{i}{pi},fibre.gamma(CS_bcc,specimenSymmetry('1'),'full'),10*degree);
            fprintf(' %s | Volume fractions | alpha fibre: %f6, gamma fibre: %f6\n',name_odf,vol_alp,vol_gam)    

        else
            plotPDF(odf{i}{pi},h{pi}(1:3),'contourf',0:0.5:4)

            disp('beta fibre')
            % look at fibre volume fractions
            vol_bet = 100 * volume(odf{i}{pi},betaFibre,10*degree);
            fprintf(' %s | Volume fractions | beta fibre: %f6\n',name_odf,vol_bet) 

        end

        mtexColorbar('title',odf_in.name{pi})
        fname = [parfile_dir '/' sample '/Figures/' name_odf prefix '_recalcPF.jpg'];
        saveFigure(fname)
        close
    end
       
end

    
    

