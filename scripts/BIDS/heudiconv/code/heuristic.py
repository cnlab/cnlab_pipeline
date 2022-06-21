import os


def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes


def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where

    allowed template fields - follow python string module:

    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    
    Example:
    t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
    func_image=create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-image_run-{item:01d}_bold')
    
    info = { t1w: [], func_image: [] }
    
        for s in seqinfo:
            if (s.dim1 == 256) and ('MPRAGE' in s.series_id):
                info[t1w].append(s.series_id)

            if (s.dim4 == 824 or s.dim4 == 813) and ('image' in s.series_id):
                info[func_image].append(s.series_id)
        return info
            
    """

    t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
    func_image=create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-image_run-{item:01d}_bold')
    
    info = {
            t1w: [], 
            func_image: []
           }
    

    for s in seqinfo:
        """
        The namedtuple `s` contains the following fields:

        * total_files_till_now
        * example_dcm_file
        * series_id
        * dcm_dir_name
        * unspecified2
        * unspecified3
        * dim1
        * dim2
        * dim3
        * dim4
        * TR
        * TE
        * protocol_name
        * is_motion_corrected
        * is_derived
        * patient_id
        * study_description
        * referring_physician_name
        * series_description
        * image_type
        """
        
            if (s.dim1 == 256) and ('MPRAGE' in s.series_id):
                info[t1w].append(s.series_id)

            if (s.dim4 == 824 or s.dim4 == 813) and ('image' in s.series_id):
                info[func_image].append(s.series_id)
            
    return info
