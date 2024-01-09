function recon_QSM_frm_MPRAGE(data_series, mask, mask_csf, B0, TE, B0_direction, max_rep, mpragepath, SAVEPATH2)
% QSM recon from MPRAGE phase, requires:
% 1) ROMEO for phase unwrapping: https://github.com/korbinian90/ROMEO
% 2) STI Suite v3.0 for background removal: https://people.eecs.berkeley.edu/~chunlei.liu/software.html
% 3) Hongfu QSM recon: https://github.com/sunhongfu/QSM/tree/master/dipole_inversion

gyro = 2*pi*42.58;       % rad/T 10^6, gyromagnetic ratio

b0_str = num2str(B0);
b0_str = strrep(b0_str,'.','p');

TE_str = num2str(TE);
TE_str = strrep(TE_str,'.','p');


for trail = 1: max_rep

    %% load MPRAGE Mag and phase
     mprage_mag_fn = [mpragepath,'mprage_mag_s',num2str(data_series),'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];
     mprage_ph_fn =  [mpragepath, 'mprage_ph_s',num2str(data_series),'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];

     nii = load_untouch_nii(mprage_mag_fn);
     mag = (single(nii.img));

     vox = nii.hdr.dime.pixdim(2:4);

     nii = load_untouch_nii(mprage_ph_fn);
     ph = (single(nii.img)); 

      %% phase unwrap using romeo

      unph_fn = [SAVEPATH2,'unph_s',num2str(data_series),'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];
      
    if ~exist(unph_fn,'file')
        parameters.output_dir = fullfile(['romeo_tmp']); % if not set pwd() is used
        parameters.mag = mag;
        parameters.mask = mask;
        parameters.calculate_B0 = false;
        parameters.phase_offset_correction = 'off';
        parameters.voxel_size = vox;
        parameters.additional_flags = '';
        mkdir(parameters.output_dir);

        [unphase, ~] = ROMEO(ph, parameters); % run
                        
        nii = make_nii(unphase,vox); save_nii(nii,unph_fn); % save

    else % load existing
        disp('unwrap file found, loaded ...')
        nii = load_untouch_nii(unph_fn); unphase = single(nii.img);
  

    end


    %% background removal, using v-sharp from STI Suite v3.0

        lfs_fn = [SAVEPATH2,'lfs_s',num2str(data_series),'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];

       [TissuePhase,NewMask] = V_SHARP(unphase, mask,'voxelsize', vox,'smvsize', 12);
       NewMask = PolishMask(NewMask);

       nii = make_nii(TissuePhase.*NewMask,vox); save_nii(nii,lfs_fn); % save

       TissuePhase = TissuePhase.*NewMask/(gyro*B0*TE*1e-3);
                            
   %% dipole inversion, using a modified version of Hongfu QSM recon to utilize GPU

       scl_factor = (7/B0) *4.44/TE; % adjust regularization based on TE & B0
       tv_reg = scl_factor*11e-4;

       disp(['data: ',num2str(data_series),' , B0: ',num2str(B0), ' TE:', num2str(TE),  ' lamda:',num2str(tv_reg)])
                    
       sus = tvdi_gpu(TissuePhase.*mask_csf,NewMask.*mask_csf, vox, tv_reg, mag, B0_direction, 500); 
         
       sus = sus.*NewMask.*mask_csf;       

       sus_fn = [SAVEPATH2,'sus_s',num2str(data_series),'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];
                    
       nii = make_nii(sus,vox); save_nii(nii,sus_fn); % save
end