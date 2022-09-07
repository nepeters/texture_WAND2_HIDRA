%% Load in MAUD ODFs
% using import tool written by Dan Savage
% https://github.com/djm87/MAUD-TBP

% main directory

parfiles     = ["Calcite_VULCAN_EWIMV10.par"];

CS = crystalSymmetry('-3m1', [4.98 4.98 17.09], 'X||a', 'Y||b*', 'Z||c*');

h = { ...
  Miller(0,1,2,CS),...
  Miller(1,0,4,CS),...
  Miller(0,0,6,CS),...
  };

for i=1:length(parfiles)
    
    parfile     = parfiles(i);
    name        = erase(parfile,["\",".par"]);
    parfile_loc = parfiles(i);
    
    % disp(['loading ODF(s) from' sample]);
    odf_in = ExtractODFFromPar(parfile_loc);

    % store in cell
    odf = {};
    % loop over phases
    for pi=1:odf_in.num

        maud_odf = odf_in.odf{pi};

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

        odf_out = strcat(name,"_mtex.odf");
        fout = fopen(odf_out,'w+');
        fprintf(fout,'alpha\tbeta\tgamma\tweight\n');
        fprintf(fout,'%g\t %g\t %g\t %g\n',maud_odf');

        % specify kernel
        psi = deLaValleePoussinKernel('halfwidth',10*degree);

        % load the ODF into the variable odf
        odf{pi} = ODF.load(odf_out,odf_in.cs{pi},odf_in.ss{pi},'density','kernel',psi,'resolution',10*degree,...
          'interface','generic',...
          'ColumnNames', { 'alpha' 'beta' 'gamma' 'weight'}, 'Matthies', 'Degrees');

        odf{pi} = rotate(odf{pi},rotation('axis',zvector,'angle',90*degree));
      
        maud_ori = orientation.byEuler(maud_odf(:,1)*degree,maud_odf(:,2)*degree,maud_odf(:,3)*degree,'ZYZ',odf_in.cs{pi},odf_in.ss{pi});

        values=odf{pi}.eval(maud_ori);
        err=sum(abs(orig_wgts-values));
        [max_err,ind]=max(abs(orig_wgts-values));
        mean_err=mean(abs(orig_wgts-values));
        fprintf(' %s | Import err : %f6, err max: %f6, err mean: %f6\n',name,err,max_err,mean_err)      

        figure
        plotPDF(odf{pi},h,'contourf',0:0.1:2.0)
        mtexColorbar('title','Calcite')
        saveFigure(char(strcat(name,"_pfs.jpg")))
        close

        t_id = textureindex(odf{pi});
        fprintf(' %s | Texture index: %f6\n',name,t_id)
        
    end
    
end

