fprintf(1,'Executing %s at %s:\n',mfilename(),datestr(now));
ver,
try,
addpath('/data00/tools/spm12mega');
con_index = 1;
cluster_forming_thr = 4.560000;
stat_filename = '/data00/projects/megameta/stanford_ks2/derivatives/nipype/task-fund_model-funded/l2analysis/stim_1.0stim_0.0/spmT_0001.nii';
extent_threshold = 0;
load '/data00/projects/megameta/stanford_ks2/derivatives/nipype/task-fund_model-funded/l2analysis/stim_1.0stim_0.0/SPM.mat'

FWHM  = SPM.xVol.FWHM;
df = [SPM.xCon(con_index).eidf SPM.xX.erdf];
STAT = SPM.xCon(con_index).STAT;
R = SPM.xVol.R;
S = SPM.xVol.S;
n = 1;

voxelwise_P_Bonf = spm_P_Bonf(cluster_forming_thr,df,STAT,S,n)
voxelwise_P_RF = spm_P_RF(1,0,cluster_forming_thr,df,STAT,R,n)

stat_map_vol = spm_vol(stat_filename);
[stat_map_data, stat_map_XYZmm] = spm_read_vols(stat_map_vol);

Z = stat_map_data(:);
Zum = Z;

        switch STAT
            case 'Z'
                VPs = (1-spm_Ncdf(Zum)).^n;
                voxelwise_P_uncor = (1-spm_Ncdf(cluster_forming_thr)).^n
            case 'T'
                VPs = (1 - spm_Tcdf(Zum,df(2))).^n;
                voxelwise_P_uncor = (1 - spm_Tcdf(cluster_forming_thr,df(2))).^n
            case 'X'
                VPs = (1-spm_Xcdf(Zum,df(2))).^n;
                voxelwise_P_uncor = (1-spm_Xcdf(cluster_forming_thr,df(2))).^n
            case 'F'
                VPs = (1 - spm_Fcdf(Zum,df)).^n;
                voxelwise_P_uncor = (1 - spm_Fcdf(cluster_forming_thr,df)).^n
        end
        VPs = sort(VPs);

voxelwise_P_FDR = spm_P_FDR(cluster_forming_thr,df,STAT,n,VPs)

V2R        = 1/prod(FWHM(stat_map_vol.dim > 1));

clusterwise_P_RF = spm_P_RF(1,extent_threshold*V2R,cluster_forming_thr,df,STAT,R,n)

[x,y,z] = ind2sub(size(stat_map_data),(1:numel(stat_map_data))');
XYZ = cat(1, x', y', z');

[u, CPs, ue] = spm_uc_clusterFDR(0.05,df,STAT,R,n,Z,XYZ,V2R,cluster_forming_thr);

clusterwise_P_FDR = spm_P_clusterFDR(extent_threshold*V2R,df,STAT,R,n,cluster_forming_thr,CPs')

,catch ME,
fprintf(2,'MATLAB code threw an exception:\n');
fprintf(2,'%s\n',ME.message);
if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
end;