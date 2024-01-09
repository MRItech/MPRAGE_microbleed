function CMB_MPRAGE_sim_func(data_series, phase, PD, sus, gm_m, wm_m, csf_m, B0, TE, vox, max_rep, SAVEPATH)
% this code uses GPU, remove 'gpuArray' if no GPU is available

%% const
ref_SNR = 32/3; % meausred at 3T, normalized by 3

siz = size(phase);

TE_str = num2str(TE);
TE_str = strrep(TE_str,'.','p');
%% initial complex MR signal

Mz_0 = gpuArray(PD.*exp(1i.*phase)) ;


%% load MRAGE recovery curves and calculate R2* estimate
switch B0

    case 1.5
            Mzr_fn = ['Mz_1p5T.mat'];
            b0_str = '1p5';
             R2s = (abs(sus).* 127 + 23)/2; % approx R2* based on suscep
    case 3
            Mzr_fn = ['Mz_3T.mat'];
             b0_str = '3';
             R2s = abs(sus).* 127 + 23;  % approx R2* based on suscep
    case 7
        Mzr_fn = ['Mz_7T.mat'];
         b0_str = '7';
             R2s = 2.28*(abs(sus).* 127 + 23) -4.23; % approx R2* based on suscep

end

% To reduce run time, load pre-computed magnetization recovery for MPRAGE readouts
% can be re=computed using: https://github.com/MRItech/mprage/blob/master/MPRAGE_phase_simulation.m

load (Mzr_fn, 'Mzr','theta') % load Mz recovery curves

Mzr = gpuArray(permute(Mzr, [2 3 1])); % move to GPU for faster computation

% masks for main tissue components 
gm_m = gpuArray(gm_m);
wm_m = gpuArray(wm_m);
csf_m = gpuArray(csf_m);

%% compute T2* and flip angle contribution factor
t2s_factor =  gpuArray(sind(theta) .* exp(-TE*1e-3.*R2s));
ks_rr = gpuArray(zeros(siz ));

N = size(phase,1); % # of slices, assuming sagital acqusition

for i=1:N   % for each k-space line along slice encoding 

    % Prepare MR signal in space domain, 
    % each tissue has its own recovery (i.e., Mzr)
    sig_r = Mz_0.*(gm_m.*Mzr(1,i) + wm_m.*Mzr(2,i)+csf_m.*Mzr(3,i) ) .*t2s_factor;  
     
    % transform into k-space
    ks = fftshift(fftn(sig_r));
    
    % store one slice encoding line in each excitation
    ks_rr(i,:,:) = ks (i,:,:); 
end

%% calcualte SNRm (magnitude SNR)

SNRm = ref_SNR*B0; % assuming linear change with B0

%% Pepare noise
Es = sum(ks_rr(:).*conj(ks_rr(:))); % MR singal energy

for trail = 1: max_rep %  repitiations with different noise instances
 
    disp(['                       trail:', num2str(trail)])

    % random complex noise
    noise = (randn(siz)) + 1j*(randn(siz));
    
    % noise energy
    En = sum(noise(:).*conj(noise(:)));

    gain = sqrt(Es/SNRm/En);
    
    % add noise in k-space 
    ks_r = (ks_rr )+ gain*(noise);

    % cut k-space corners, i.e., elipitcal filtering
    ks_r = elept_filt(ks_r);

    % transform into space domain
    sign_rcv = gather(ifftn(ifftshift(ks_r)));

% store MPRAGE magnitude and phase 
MPRAGE_mag = single(abs(sign_rcv));
MPRAGE_phs = single(angle(sign_rcv));


mprage_mag_fn = [SAVEPATH,'mprage_mag_s',data_series,'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];
mprage_ph_fn =  [SAVEPATH, 'mprage_ph_s',data_series,'_itr',num2str(trail),'_b',b0_str,'_TE',TE_str,'.nii'];

nii = make_nii(MPRAGE_mag,vox); save_nii(nii,mprage_mag_fn);
nii = make_nii(MPRAGE_phs,vox); save_nii(nii,mprage_ph_fn);

clear  ks_r  noise sign_rcv MPRAGE_mag MPRAGE_phs nii
end



end

% k-space eliptical filter
function y = elept_filt(x)
N = size(x);
[ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);
a = N(1)/2;
b = N(2)/2;
c = N(3)/2;

mask = zeros(N);
mask((kx/a).^2 + (ky/b).^2 + (kz/c).^2 <=1) = 1;
y = x.*mask;
end
