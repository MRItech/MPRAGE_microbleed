function measure_suscep_and_volume(data_series, sus_org, meas_mask, B0,TE, max_rep, SAVEPATH2)
% measure microbleed size and susceptibility


b0_str = num2str(B0);
b0_str = strrep(b0_str,'.','p');

TE_str = num2str(TE);
TE_str = strrep(TE_str,'.','p');

%% load data

for trail = 1: max_rep

    sus_fn = [SAVEPATH2,'sus_s',num2str(data_series),'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];
    nii = load_untouch_nii(sus_fn); sus = single(nii.img);
    siz = size(sus);

    vox = nii.hdr.dime.pixdim(2:4);
    vox_vol = prod(vox);            % voxel volume
    vox_area = prod(vox(1:2));      % voxel area
     %% arrays to store measurements

    max_sus = zeros(2,20,'single'); % max sus within ROI
    vol = max_sus;                  % ROI volume
    area = max_sus;                 % ROI area
    msus_vol = max_sus;             % mean suscep over volume
    msus_area= max_sus;             % mean suscep over area
    tsus_vol = max_sus;             % total suscep within volume
    tsus_area =max_sus;             % total suscep over area
   
    roi_mask_all = zeros(siz,'single');

    % loop for all microbleeds
    for msk =1: 20
        
        roi_mask = zeros(siz,'single');
        roi_mask(meas_mask==msk) = 1;
        
        % sus mean & max
        sus_roi = sus.*roi_mask;
        max_sus(1,msk) = max(sus_roi(:)); % max suscep within mask

        org_sus_roi = sus_org.*roi_mask;
        max_sus(2,msk) = max(org_sus_roi(:)); % max in ground truth
        
        % keep voxels > 0.2*max && > 0.05
        HM_sus = sus_roi; 
        HM_sus(HM_sus<max_sus(1,msk)/5) = 0;
        HM_sus(HM_sus<0.05) = 0;
        HM_msk = HM_sus;  HM_msk(HM_msk>0) = 1;


        HM_org_sus = org_sus_roi; 
        HM_org_sus(HM_org_sus<max_sus(2,msk)/5) = 0;
        HM_org_sus(HM_org_sus<0.05) = 0;
        HM_org_msk = HM_org_sus;  HM_org_msk(HM_org_msk>0) = 1;
                    
        % store the produced ROI mask
        roi_mask_all = roi_mask_all + msk*HM_msk; 

         % volume  at 20% max
         vol(1,msk) =sum(HM_msk(:)) ;
         vol(2,msk) =sum(HM_org_msk(:)) ;
         
         % find middle slice and calculate its area
         ar = sum(sum(HM_msk,1),2) ;
         [area(1,msk), indx] = max(ar);
         msk_area = zeros(siz,'single');
         msk_area(:,:,indx) =  HM_msk(:,:,indx);

         ar = sum(sum(HM_org_msk,1),2) ;
         [area(2,msk), indx] = max(ar);
         msk_org_area = zeros(siz,'single');
         msk_org_area(:,:,indx) =  HM_org_msk(:,:,indx);

         % mean suscep over volume and area
         msus_vol(1,msk) = sum(HM_msk(:).*sus_roi(:))/sum(HM_msk(:));
         msus_area(1,msk) = sum(msk_area(:).*sus_roi(:))/sum(msk_area(:));
                    
         msus_vol(2,msk) = sum(HM_org_msk(:).*org_sus_roi(:))/sum(HM_org_msk(:));
         msus_area(2,msk) = sum(msk_org_area(:).*org_sus_roi(:))/sum(msk_org_area(:));
        
         % total suscep within volume and area

         vol= vol *vox_vol; 
         area = area*vox_area;

        tsus_vol(1,msk) = msus_vol(1,msk)*vol(1,msk);
        tsus_area(1,msk)= msus_area(1,msk)*area(1,msk);

        tsus_vol(2,msk) = msus_vol(2,msk)*vol(2,msk);
        tsus_area(2,msk)= msus_area(2,msk)*area(2,msk);

    end

    %% store results
         mat_fn = [SAVEPATH2,'meas_s',num2str(data_series),'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.mat'];
         save(mat_fn,'max_sus','vol','area','msus_vol','msus_area','tsus_area','tsus_vol')
                    
         msk_fn = [SAVEPATH2,'measROI_s',num2str(data_series),'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];
         nii = make_nii(roi_mask_all); save_nii(nii,msk_fn)


end