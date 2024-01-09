% Simulation code of quantifying microbleeds using MPRAGE-QSM
% Paper: Naji N, et al. Quantifying cerebral microbleeds using MPRAGE-QSM,
% NMR in Biomed

clear
%% General settings
update_MBs_locations = 0;    % demands huge memory of min. 64GB
recompute_MPRAGE = 0;        % recompute MPRAGE complex data
update_recon = 0;            % re-caclulate QSM
do_measurements = 0;         % measure MB volume and susceptibility 

%% Paramters to set
d_sus = [0.05 .1 .2 .4 .7 1 1.5 2];      % ppm, susceptibility levels
d_diameter = [ 0.3 0.5 0.9 1.7 3.5 6.9]; % mm, microbleeds dimeters
B0 = [ 1.5 3 7]; % T, field strength

TE = [ 2.37, 4.44]; %[ 2.37, 3 , 3.5, 4.44]; % ms, echo-times

max_rep = 10; % # of experiment repitiations with different random noise

B0_direction = [ 0 0 1]; % B0 direction
vox = [0.1 0.1 0.1;];    % mm, new voxel size
gyro = 2*pi*42.58;       % rad/T 10^6, gyromagnetic ratio

datapath_src_low = ['low_res_rsc/']; 
SAVEPATH = 'mprage_out_lp/';
SAVEPATH2 = 'recon_out_lp/';
sus_fn = ['org/','head_sus.nii'];

%% laod original data

hi_res_sus_file = ['Hi_res/sus_inp.nii'];

if ~exist(hi_res_sus_file,'file') % interpolate only once
    
    disp('======= Base susceptibility dist. of 0.1mm res was not found, generating ... ')
    % suscep distirbution
    nii = load_untouch_nii(sus_fn);
    head_sus = single(nii.img);

    imsiz_org = size(head_sus);

    % original voxel size
    vox_org = nii.hdr.dime.pixdim(2:4);
    %% up sample to 0.1 mm

    % by x3.3 factor, from 0.33 mm to 0.1 mm
    sus_intp = imresize3(head_sus,3.3);

    % zeropadding in image space  to faciliate down-sampling by integer
    sus_intp = padarray(sus_intp,[0 1 1],0,'post'); % one side
    sus_intp = padarray(sus_intp,[52 246 160],0,'both'); % both sides

   
    imsiz = size(sus_intp); % array size

    % save
    nii = make_nii((sus_intp),vox);
    mkdir('Hi_res')
    save_nii(nii,hi_res_sus_file)

    clear sus_intp head_sus nii

else % get image size
    disp('======= Base susceptibility dist. of 0.1mm exists, skip generation setp ... ')
    hdr = load_nii_hdr(hi_res_sus_file);
    imsiz = hdr.dime.dim(2:4);
end


%% Generate Dipole kernel at 0.1 mm resolution

hi_res_Dip_file = ['Hi_res/D_0p1res.nii'];

if ~exist(hi_res_Dip_file,'file') % generate only once, demands huge memory
     disp('======= Dipole kernel of 0.1mm res was not found, generating ... ')

     D = dipole_kernel_2(imsiz, vox, B0_direction,'kspace' );
     nii = make_nii(single(D),vox);
     save_nii(nii,hi_res_Dip_file)

     clear nii D
end
 
%% Add microbleeds to the high resolution susceptibility map & produce total field shift

coord= load ('cmb_locations.mat'); % coordinates  of MBs locations

if update_MBs_locations
    disp('======= Updating microbleeds locations and re-generating high-res field shift... ')

    mkdir(datapath_src_low)
    add_MB_at_hiRes_calc_field(d_sus,d_diameter,coord,vox,imsiz,hi_res_sus_file,hi_res_Dip_file)
end


%%  Proton density map & tissue segmentation masks

PD_fn = ['org/PD_MAP.nii'];
nii = load_untouch_nii(PD_fn);
PD = (single(nii.img));             % synthatic proton density map

siz = size(PD);
vox = nii.hdr.dime.pixdim(2:4);

nii = load_untouch_nii('org/seg_mask.nii'); 
seg_mask2 = (single(nii.img)); 


gm_m = zeros(siz,"single");
gm_m(seg_mask2==2)=1;       % GM mask


wm_m = zeros(siz,"single");  
wm_m(seg_mask2==3)=1;       % WM mask

csf_m  = zeros(siz,"single");  
csf_m(seg_mask2==1)=1;      % CSF mask

%%  == loop for MPRAGE simulation =========

[B0_list, TE_list] = meshgrid(B0, TE); % all (B0,TE) pairs for simulation

B0_list = B0_list(:);
TE_list = TE_list(:);

if recompute_MPRAGE
    disp('======= Simulating MPRAGE magntiude and phase with additive noise...')
    mkdir(SAVEPATH)
    for data_series  = 1:6 % Loop for different microbleeds distibutions

        disp(['Prcoess data:', num2str(data_series),' .....',])
        tfs_fn = [datapath_src_low,'tfs_lp',num2str(data_series),'.nii'];
        sus_fn = [datapath_src_low,'sus_r',num2str(data_series),'.nii'];

        nii = load_untouch_nii(tfs_fn);
        tfs = (single(nii.img));        % total field shift

        nii = load_untouch_nii(sus_fn); 
        sus = (single(nii.img));        % true susceptibiltiy 

        % loop through different TE and B0 

        for sett = 1: length(B0_list)

            % form phase at TE_i and B0_i
            disp(['            B0:', num2str(B0_list(sett)),'  TE: ',num2str(TE_list(sett)),' .....',])

            phase = tfs*gyro*B0_list(sett)*TE_list(sett)*1e-3; % rad
        
            % call function to simulate MPRAGE magnitude and phase 
            CMB_MPRAGE_sim_func(num2str(data_series),phase,PD,sus,gm_m,wm_m,csf_m,B0_list(sett),TE_list(sett),vox, max_rep, SAVEPATH);
        end

    end

end

%%  ============ QSM reconstruction ==================

    % brain mask
    mask_fn = ['org/mask_r.nii'];
    nii = load_untouch_nii(mask_fn);
    mask = (single(nii.img)); 

   % mask out CSF
    mask_csf = 1-csf_m; 

    %% loop dataset
if  update_recon 

    disp('======= Reconstruction of QSM from MPRAGE Phase...')
    mkdir(SAVEPATH2)

    for data_series  =1:6
        % loop through different TE and B0
        for sett = 1: length(B0_list)
        
            recon_QSM_frm_MPRAGE(data_series, mask, mask_csf, B0_list(sett),TE_list(sett), B0_direction,  max_rep, SAVEPATH, SAVEPATH2)
        

        end

    end
end
%% ================ detection and measurements =================

if do_measurements
    disp('======= Measuring volume and susceptibiltiy of microbleeds...')
    for data_series  =1:6
        %         % load measurements mask
        fn_meas_mask = [datapath_src_low,'meas_mask_',num2str(data_series),'.nii'];
        nii = load_untouch_nii(fn_meas_mask);
        meas_mask = (single(nii.img)); 


        % laod true susceptibiltiy
        org_sus_mask = [datapath_src_low,'sus_r',num2str(data_series),'.nii'];
        nii = load_untouch_nii(org_sus_mask);
        sus_org = (single(nii.img)); 
        
        % loop through different TE and B0
        for sett = 1: length(B0_list)

            measure_suscep_and_volume(data_series, sus_org, meas_mask, B0_list(sett),TE_list(sett), max_rep, SAVEPATH2)
        end
    end
end



%% ============= Present Results ===================

% Fig 1: sus, mag, phase

% groundtruth 
nii = load_untouch_nii([datapath_src_low,'sus_r3.nii']);
resl_{1}.img  = permute(single(nii.img),[2 1 3]);

nii = load_untouch_nii([SAVEPATH,'mprage_ph_s3_itr1_b3_TE4p44.nii']);
resl_{2}.img  = permute(single(nii.img),[2 1 3]);

nii = load_untouch_nii([SAVEPATH,'mprage_mag_s3_itr1_b3_TE4p44.nii']);
resl_{3}.img  = permute(single(nii.img),[2 1 3]);

slices = [111, 128, 146];
tit_str = {'True Susceptibility', 'MPRAGE Phase', 'MPRAGE magnitude'};

tiledlayout(3,3,'TileSpacing','none');

gray_scle = [-.1 .15; -3.14 3.14;.006 .03];
for ro = 1:3
    for col = 1: 3
        nexttile
         imshow(resl_{ro}.img(:,:,slices(col)),gray_scle(ro,:))

         if col ==1
             ylabel(tit_str{ro})
         end

         if col ==3
             colorbar
         end

    end
end

clear resl_

%% Fig 6: recon suscep at diff B0, TE

figure
% recon at 1.5T
nii = load_untouch_nii([SAVEPATH2,'sus_s2_itr1_b1p5_TE2p37.nii']);
sus_{1}.img = permute(single(nii.img),[2 1 3]);

nii = load_untouch_nii([SAVEPATH2,'sus_s2_itr1_b1p5_TE4p44.nii']);
sus_{2}.img  = permute(single(nii.img),[2 1 3]);

% recon at 3T

nii = load_untouch_nii([SAVEPATH2,'sus_s2_itr1_b3_TE2p37.nii']);
sus_{3}.img  = permute(single(nii.img),[2 1 3]);

nii = load_untouch_nii([SAVEPATH2,'sus_s2_itr1_b3_TE4p44.nii']);
sus_{4}.img = permute(single(nii.img),[2 1 3]);

% recon at 7T
nii = load_untouch_nii([SAVEPATH2,'sus_s2_itr1_b7_TE2p37.nii']);
sus_{5}.img  = permute(single(nii.img),[2 1 3]);

nii = load_untouch_nii([SAVEPATH2,'sus_s2_itr1_b7_TE4p44.nii']);
sus_{6}.img  = permute(single(nii.img),[2 1 3]);

% groundtruth 
nii = load_untouch_nii([datapath_src_low,'sus_r2.nii']);
sus_{7}.img  = permute(single(nii.img).*mask,[2 1 3]);

% plot
slices = [48, 96, 110, 148];
tit_str = {'1.5T, TE:2.37ms', '1.5T, TE:4.44ms', '3T, TE:2.37ms', '3T, TE:4.44ms','7T, TE:2.37ms','7T, TE:4.44ms', 'True X'};

tiledlayout(4,7,'TileSpacing','none');
for ro = 1:4
    for col = 1: 7
        nexttile
        if col <3 % 1.5 T
             imshow(sus_{col}.img(:,:,slices(ro)),[-.3 .4])
        else
             imshow(sus_{col}.img(:,:,slices(ro)),[-.1 .15])
        end

        if ro ==1
            title(tit_str{col})
        elseif ro ==4
            colorbar('southoutside')
        end
        

    end
end
clear sus_