% add microbleeds to the high resolution susceptibility map & produce total field shift
function add_MB_at_hiRes_calc_field(d_sus,d_diameter,coord,vox,imsiz,hi_res_sus_file,hi_res_Dip_file)
% max of 20 mcirobleeds per susceptibilty distribution
% since this function is memory demanding, neasted functions were used to
% reduce the required memory.

[amp,d ] = meshgrid(d_sus,d_diameter);

ln_amp = length(amp(:)); % # of pairs

% make the # of pairs integer multiple of 20
ln_amp_20 = ceil(ln_amp/20)*20;

amp = [amp(:); zeros(ln_amp_20-ln_amp,1)];
d = [d(:); zeros(ln_amp_20-ln_amp,1)];

amp = [amp(:);amp(:)]; % twice, once in cortical and once in subcortical
d =[d(:);d(:)];        % twice, once in cortical and once in subcortical

full_loc = zeros(20,3,'single'); % coordinates of 20 locations


full_loc(1:2:end,:) = coord.locations_edge; % odd, cortical 
full_loc(2:2:end,:) = coord.locations_mid;  % even, subcortical

% total of 48 (sus, diameter) pairs, at least 3 sets of max 20 microbleeds
full_loc = [full_loc;full_loc;full_loc;full_loc(end:-1:1,:);full_loc(end:-1:1,:);full_loc(end:-1:1,:)];

% loop to add microbleeds
sr =1;
for l=1:20:length(amp)
    
        amp_o = amp(l:l+20-1);
        d_o = d(l:l+20-1);
        loc = full_loc(l:l+20-1,:);

        sus_intp = load_untouch_nii(hi_res_sus_file);
        sus_intp = single(sus_intp.img);
%     N = size(sus_intp);
    add_cmb_multiple
    sr = sr+1;
end



clear sus_intp 

%% ========================= functions ================
function  add_cmb_multiple

mask = zeros(size(sus_intp), 'uint8');
cmb_no = length(amp_o);

        for i = 1: cmb_no
            % form cmb with amp and d
            cmb = form_cmb(d_o(i));
            
            if ~isempty(cmb)
            % add at location by masking and addition
            sz = (size(cmb)-1)/2;
            tmp = sus_intp(loc(i,1)-sz(1):loc(i,1)+sz(1),loc(i,2)-sz(2):loc(i,2)+sz(2),loc(i,3)-sz(3):loc(i,3)+sz(3)) ;
            tmp = tmp.* (1-cmb) + cmb.*amp_o(i);
            sus_intp(loc(i,1)-sz(1):loc(i,1)+sz(1),loc(i,2)-sz(2):loc(i,2)+sz(2),loc(i,3)-sz(3):loc(i,3)+sz(3)) = tmp;
            
            sz_all(i,:) = size(cmb);
            loc_all(i,:) = loc(i,:);

            else % zero diameter
                sz_all(i,:) = [0, 0, 0]; 
                loc_all(i,:) = loc(i,:);
            end

            % mask
            if d_o(i) < 3.5
                cmb2 = form_cmb(3.5);
            else
                cmb2 = form_cmb(d_o(i)*2 -0.1);
            end

            if ~isempty(cmb2)
                sz = (size(cmb2)-1)/2;
                mask(loc(i,1)-sz(1):loc(i,1)+sz(1),loc(i,2)-sz(2):loc(i,2)+sz(2),loc(i,3)-sz(3):loc(i,3)+sz(3)) = uint8(cmb2*i);
            end

        end
        disp('CMBs added ..')
        % save mask after down-sampling
        mask_lowres = imresize3(padarray(mask,[0 0 1],0,'pre'),0.117647,'nearest'); % 0.1mm/0.85mm
        clear mask

        % SAVE modifed sus
        nii = make_nii(sus_intp(53:end-52,246+1:end-246,160+1:end-160),vox);
        save_nii(nii,  ['Hi_res/','sus_',num2str(sr),'.nii'] )
        clear nii
        % save lower resolution version
        sus1 = imgaussfilt3(sus_intp,4.25,'FilterDomain','spatial');
        disp('gauss filtered ...')
        FOV = (imsiz + [0 0 1]).*vox;
        ima_lowres = imresize3(padarray(sus_intp,[0 0 1],0,'pre'),0.117647); % 0.1mm/0.85mm
        new_size = size(ima_lowres);
        new_vox = FOV./new_size;
        nii = make_nii(ima_lowres(:,1:end-1,:),new_vox);
        save_nii(nii,['low_res_rsc/','sus_r',num2str(sr),'.nii'])
        
        nii = make_nii(mask_lowres(:,1:end-1,:),new_vox);
        save_nii(nii,['low_res_rsc/','meas_mask_',num2str(sr),'.nii'])

        ima_lowres = imresize3(padarray(sus1,[0 0 1],0,'pre'),0.117647); % 0.1mm/0.85mm
        nii = make_nii(ima_lowres(:,1:end-1,:),new_vox);
        save_nii(nii,['low_res_rsc/','sus_gr',num2str(sr),'.nii'])

        clear nii ima_lowres sus1 mask_lowres
        
        % high res field calculation in freq domain
        DC_centre                 = floor(imsiz/2)+1;
        lower_bound               = DC_centre - ceil(imsiz.*0.117647/2) ;
        upper_bound               = lower_bound      + ceil(imsiz.*0.117647) - 1; upper_bound(2) = upper_bound(2) +1;
        new_imsiz               = upper_bound-lower_bound+1; scale_f = imsiz./new_imsiz;
        smoth_width                  = 30; % # of voxels where the filter will drop off from 1 to 0.
        [kfx,kfy,kfz]           = ndgrid(m1_cosfilt(new_imsiz(1),smoth_width),m1_cosfilt(new_imsiz(2),smoth_width),m1_cosfilt(new_imsiz(3),smoth_width));
        kf                      = min(min(kfx,kfy),kfz); clear kfx kfy kfz
        % simulate field shift
            sus_intp = fftn(sus_intp);
            % split real and imag to facililitate memory usage
            sus_r = real(sus_intp);        
            sus_i = imag(sus_intp);
             sus_intp = [];
             disp('fft done..')

             % load calculated diple kernel
             D = load_untouch_nii(hi_res_Dip_file);
             D = single(D.img);

             sus_r = sus_r.*D;
             sus_i = sus_i.*D;
             clear D
             disp('dipole filtering done ..')

             
             sus_r = sus_r + 1j*sus_i;
             sus_i = [];

             sus_i = fftshift(sus_r);
             sus_i = sus_i(lower_bound(1):upper_bound(1),lower_bound(2):upper_bound(2),lower_bound(3):upper_bound(3)).*kf;

             sus_i = ifftshift(sus_i);
             ima_lowres2 = real(ifftn(sus_i))./prod(scale_f);
             
             % save low-res field shift
            nii = make_nii(ima_lowres2(:,2:end-1,:),new_vox); clear ima_lowres2
            save_nii(nii,['low_res_rsc/','tfs_lp',num2str(sr),'.nii']); clear nii
            
            sus_r = real(ifftn(sus_r));
            tfs = sus_r;
            sus_r = [];
            disp('ifft done ..')
        

            clear nii ima_lowres 


        % crop hi-res field to reduce file size 
        tfs = tfs(53:end-52,246+1:end-246,160+1:end-160);
       
        % save hi-res field shift
        nii = make_nii(tfs,vox);
        save_nii(nii,  ['Hi_res/','tfs_',num2str(sr),'.nii'] )
        clear nii tfs
        save(['Hi_res/','series_',num2str(sr),'.mat'],"loc_all","sz_all",'amp_o','d_o')

   
end

function cmb = form_cmb(d)
% generate shperical microbleeds with specificed diameter
cmb = [];
switch single(d)
    case 0
        return;  % zero diameter
    case single(0.3)
        tmp = strel('sphere',1);
    case single(0.5)
        tmp = strel('sphere',2);   
    case single(0.9)
        tmp = strel('sphere',4);
    case single(1.7)
        tmp = strel('sphere',8);
    case single(3.5)
        tmp = strel('sphere',17);
    case single(6.9)
       tmp = strel('sphere',34); 
    case single(13.7)
       tmp = strel('sphere',68);
end

cmb = single(tmp.Neighborhood);
end
end