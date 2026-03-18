import os, os.path as op
import sys, pandas as pd
import re, glob, argparse, subprocess, shutil
import datalad.api as dl, json

"""

Usage


    Pre CuBIDS:

    python data-bidsify.py \
        <path-to-working-directory> \
        -d /path/to/dicoms \
        -dd list_of_subs_dicom_directories \
        -i /path/to/store_dicoms \
        -o /path/to/BIDS \
        -s subject \
        -ss list_of_sessions

    Post CuBIDS:

    python data-bidsify.py \
        <path-to-working-directory> \
        -o /path/to/BIDS \
        -s subject \
        -ss list_of_sessions \
        -m /path/to/cubids_apply_script.sh

"""

bgd = '/mnt/isilon/bgdlab_processing/'
simg = '/mnt/isilon/bgdlab_processing/code/containers/clin-dataorg-0.0.2.sif'

def print_cmd(result):

    '''
    Print the command used for a subprocess run

    Args:
        result: Output of subprocess run
    '''
    fstr =" ".join(result.args).replace("'","")
    print(fstr)


def is_clean(path):
    
    '''
    Check if datalad dataset is clean.
    '''

    result = dl.status(path = path, dataset = path)

    for file_info in result:
        print("state", file_info['state'])
        if file_info['state'] == 'clean':
            return True
        else:
            return False
     

def create_config(path):

    '''
    Create dcm2niix config file used to run heudiconv.
    
    Args:
        path: Path to save config file.
    
    Return:
        path to new config file
    '''

    dcmconfig = op.join(path, 'dcm2niix_config.json')
    if not op.exists(dcmconfig):
        with open(dcmconfig, 'a') as f:
            f.write('{"m": "y"}')   
    return dcmconfig       


def readJson(path, key):

    '''
    Read a json file to extract a key value
    
    Args:
        path: Full path to json file
        key: Key in json file to lookup its value

    Return:
        Value belonging to input key or None if key not found
    '''

    lst = []
    with open(path) as f:
        f = json.load(f)
        if key in list(f):
            lst.append(f[key])

    if len(set(lst)) == 1:
        return lst[0]
    elif len(set(lst)) > 1:
        return ",".join(filter(bool, lst))
    else:
        return None


def verify_output(subj, sess, niftis, jsons):
    
    '''
    Check output of heudiconv processing for unpaired nifti/jsons, 
    non-brain scans, errant DWI scans, and multi-echo scans.

    Args:
        subj: Subject ID
        sess: Session ID
        niftis: All niftis in a subject's session directory 
        jsons: All jsons in a subject's session directory
    
    Return:
        remove_files: List of files to be deleted
        change_files: List files to be renamed or moved to problematic_scans directory
    '''

    remove_files = []
    change_files = []

    for n in niftis:
        basename = n.split(".nii")[0]
        pairedJson = [i for i in jsons if (basename in i)]
        prefix = f"sub-{subj}_ses-{sess}"

        # Remove nifti with no paired json
        if not pairedJson:

            remove_files.append([n,"remove"])

        else:

            pairedJson = pairedJson[0]
            echoNum = None if not readJson(pairedJson, "EchoNumber") else readJson(pairedJson, "EchoNumber")
            descrip = 'n/a' if not readJson(pairedJson, "SeriesDescription") else readJson(pairedJson, "SeriesDescription").lower()
            protocol = 'n/a' if not readJson(pairedJson,"ProtocolName") else readJson(pairedJson,"ProtocolName").lower()
            stepDescrip = 'n/a' if not readJson(pairedJson, "ProcedureStepDescription") else readJson(pairedJson, "ProcedureStepDescription").lower()
            bodypart =  'n/a' if not readJson(pairedJson, "BodyPartExamined") else readJson(pairedJson, "BodyPartExamined").lower()
    

            # Remove non-brain scans
            exclude = ["heart","chest",'temporomandibular']
            res = [i for i in exclude if(i in bodypart or i in stepDescrip or i in descrip or i in protocol)]
            if res:
                print(f"IDENTIFIED A NON BRAIN SCAN: {res}")
                
                n2 = n.replace("BIDS","problematic_scans")
                change_files.append([n,n2])
                j2 = pairedJson.replace("BIDS","problematic_scans")
                change_files.append([pairedJson,j2])
                print(f"moved {n} to {n2}")
            
            else:
            
                # Check dwi data 
                if '/dwi' in n:
                    _bvec = n.replace("nii.gz","bvec")
                    _bval =n.replace("nii.gz","bval")
                    if not op.exists(_bval) or not op.exists(_bvec):
                        remove_files.append([n,"remove"])
                        remove_files.append([pairedJson,"remove"])
                        if op.exists(_bvec):
                            remove_files.append([_bvec,"remove"])
                        else:
                            remove_files.append([_bval,"remove"])
                
                else:
                    # Check for multi-echo data
                    if re.findall(r'w[0-9]+.|heudiconv', n):

            
                        if not echoNum:

                            # For scans with duplicate series numbers, segregate them from the rest of the data
                            n2 = n.replace("BIDS","problematic_scans")
                            change_files.append([n,n2])
                            j2 = pairedJson.replace("BIDS","problematic_scans")
                            change_files.append([pairedJson,j2])
                            print(f"moved {n} to {n2}")
                            

                        else:

                            print(f"\nechoNum: {echoNum}")
                            # For multi-echo data, add echo to file name
                            echo = "echo-" + str(echoNum)      
                            n2 = generate_bids_key(n, prefix, echo)
                            n2 = re.sub('[0-9]+.nii.gz','.nii.gz', n2)
                            print(f"moved {n} to {n2}")
                            change_files.append([n,n2])
                            j2 = n2.replace("nii.gz","json")
                            change_files.append([pairedJson,j2])
                            

                    else:
                        if "dual_echo" in descrip or "dual_echo" in protocol or "pd+" in descrip or "pd+" in protocol:

                            echo = "echo-2" if not echoNum else ("echo-" + str(echoNum))  
                            n2 = generate_bids_key(n, prefix, echo)   
                            n2 = re.sub('[0-9]+.nii.gz','.nii.gz',n2)
                            print(f"moved {n} to {n2}")
                            change_files.append([n,n2])
                            j2 = n2.replace("nii.gz","json")
                            change_files.append([pairedJson,j2])
                        
    for j in jsons:

        # Remove json with no paired nifti
        basename = j.split(".json")[0]
        pairedNifti = [i for i in niftis if (basename in i)]
        if not pairedNifti:
            remove_files.append([j,"remove"])

    return remove_files, change_files


def get_bids_key_value(key, fn):

    '''
    Lookup a key-value pair in a BIDS filename

    Args:
        key: BIDS key to lookup value
        fn: Bidsified nifti file 
    Return:
        value: value belonging to BIDS key
    '''
    
    value = str(fn.split(key)[1].split("_")[0])
    return value


def generate_bids_key(fn, prefix, toInsert):
    
    '''
    Insert value into a filename maintaining BIDS complaince 

    Args:
        fn: Nifti filename
        toInsert: value to insert into nifti filename
        prefix: sub-<>_ses-<>

    Return:
        wholename: Updated filename with inserted key-value pair
    '''

    suffix_parts = [
        None if (not "sample-" in fn and not "sample-" in toInsert) else "sample-fetal",
        None if (not "acq-" in fn and not "acq-" in toInsert) else toInsert if "acq-" in toInsert else "acq-%s" % get_bids_key_value("acq-",fn),
        None if (not "ce-" in fn and not "ce-" in toInsert) else "ce-gad",
        None if (not "dir-" in fn and not "dir-" in toInsert) else toInsert if "dir-" in toInsert else "dir-%s" % get_bids_key_value("dir-",fn),
        None if (not "rec-" in fn and not "rec-" in toInsert) else toInsert if "rec-" in toInsert else "rec-%s" % get_bids_key_value("rec-",fn),
        None if (not "inv-" in fn and not "inv-" in toInsert) else toInsert if "inv-" in toInsert else "inv-%s" % get_bids_key_value("inv-",fn),
        None if (not "run-" in fn and not "run-" in toInsert) else toInsert if "run-" in toInsert else "run-%s" % get_bids_key_value("run-",fn),
        None if (not "echo-" in fn and not "echo-" in toInsert) else toInsert if "echo-" in toInsert else "echo-%s" % get_bids_key_value("echo-",fn),
        None if (not "mt-" in fn and not "mt-" in toInsert) else "mt-on",
        None if (not "part-" in fn and not "part-" in toInsert) else toInsert if "part-" in toInsert else "part-%s" % get_bids_key_value("part-",fn),
        re.findall(r'T1w|T2w|T2starw|dwi|PDT2|ep2d|FLAIR|PDw|unknown',fn)[0]
    ]

    # filter those which are None, and join with _
    suffix = "_".join(filter(bool, suffix_parts))
    newname = prefix + '_' + suffix
    basename = op.basename(fn).split('.nii')[0]
    outFn = fn.replace(basename, newname)
    
    return outFn


def correction(subj, sess, lst, df=False):

    '''
    Delete, rename, or move scans from the BIDS directory

    Args:
        subj: Subject name
        sess: Session name
        lst: List of files to either delete or rename or move to another directory
        df: Dataframe to keep track of scans moved to the problematic directory

    Return:
        made_changes: True/False if files were deleted, renamed, or moved
    '''

    subj = 'sub-' + subj
    sess = "ses-" + sess
    made_changes = False
    
    for group in lst:

        a = group[0]
        b = group[1]
        if b == "remove":
            if op.exists(a):
                util = f'rm {a}'
                print(util)
                os.system(util)
                made_changes = True
        else:
            if 'heudiconv' in b and not 'problematic' in b:
                b = b.split("_heudiconv")[0] + re.findall(r'.nii.gz|.json|.bvec|.bval', b)[0]
            basepath = op.dirname(b)

            if op.exists(a) and op.exists(basepath) and not op.isfile(b):

                util = f"mv {a} {b}"
                made_changes = True
                
            elif op.isfile(b):

                print(f"{b} already exists")
                suffix = re.findall(r'nii.gz|json|bvec|bval', a)[0]
                files = [i for i in os.listdir(op.dirname(a)) if (suffix in i)]
                b = iterate_run(b, files)
                util = f"mv {a} {b}"
                made_changes = True

            elif not op.exists(basepath):

                print(f"{basepath} directory doesn't exist. Creating it . . .")
                os.system(f"mkdir -p {basepath}")
                util = f"mv {a} {b}"
                made_changes = True

            elif not op.exists(a):
                pass # TO DO

        print(util)
        os.system(util)

        if "problematic" in b:
            df.loc[len(df)] = [subj, sess, b]

    return made_changes


def check_output(wrkDir, subj, sess, dicomSets, outSesDir, probs, outDf, relaunch=False):
    
    '''
    Compare heudiconv processing logs and BIDS subject session directory
    to verify that dicoms passing inclusion criteria become nifti(s).
    And updates a dataframe to keep track of heudiconv's performance

    Args:
        subj: Subject name
        sess: Session name
        dicomSets: 1 level down from the GCP directory, equal to the study UID
        outSesDir: output subject session directory (e.g. /BIDS/sub-001/ses-a)
        probs: Dataframe to keep track of scans moved to the problematic directory
        outDf: Dataframe to keep track of heudiconv processing performance 
        relaunch: True/False if heudiconv was previously relaunched

    Return:
        True/False if heudiconv performed successfully
    '''
    
    resubmit = False
    expected_series_n = [i.rsplit("/",1)[0] for i in glob.glob(f"{dicomSets}/**/*.dcm")]
    target_n = len(list(set(expected_series_n)))
    # target_n = len(glob.glob(f"{dicomSets}/*"))
    
    # Get information on how heudiconv performed
    fstr = f"/logs/heudiconv/bidsify_{subj}_{sess}.txt" if relaunch is False else f"/logs/heudiconv/bidsify_redo_{subj}_{sess}.txt"
    
    print(relaunch)
    print(wrkDir,fstr)
    slurm_logs = glob.glob(wrkDir+fstr)
    print(slurm_logs)
    

    if relaunch is False:
        actual_series_n, expected_scans_n, anat_n, dwi_n, performance = fetch_answers(slurm_logs, target_n) 
    else:
        actual_series_n, expected_scans_n, anat_n, dwi_n, performance = fetch_answers(slurm_logs, target_n, True) 

    '''
    If the BIDS directory was not created, it could be because:
    1. an error with heudiconv caused it to stop or
    2. no data passed through the heuristic filters
    If reason #1, submit for series wise heudiconv runs
    If reason #2, resubmit if heudiconv missed registering all series, if not don't resubmit
    '''
    if not op.exists(outSesDir):
        if performance == 'fail':
            print("Heudiconv failed.")
            if relaunch is False:
                print("\tResubmitting . . .")
                outDf.loc[len(outDf)] = create_dict(subj, sess, target_n, actual_series_n, 0, 0, 0, 0, performance)
                return True, False
            else:
                outDf.loc[(outDf.subject_id == f"sub-{subj}")&(outDf.session_id == f"ses-{sess}"),'final_run'] = performance
                return False, False
        else:
            print("No data passed through heuristic filter.")
            if relaunch is False:
                if actual_series_n < target_n:
                    resubmit = True
                    performance = 'missed'
                    print("\tResubmitting . . .")
                else:
                    performance = 'success'
                outDf.loc[len(outDf)] = create_dict(subj, sess, target_n, actual_series_n, 0, 0, 0, 0, performance)
            else:
                if actual_series_n < target_n:
                    resubmit = True
                    performance = 'missed'
                else:
                    performance = 'success'
                outDf.loc[(outDf.subject_id == f"sub-{subj}")&(outDf.session_id == f"ses-{sess}"),'final_run'] = performance
            return resubmit, False
    else:
        os.system(f"chmod +rw -R {outSesDir}")
    
    actual_anat_n = len(glob.glob(outSesDir+"/anat/*nii.gz"))
    actual_dwi_n = len(glob.glob(outSesDir+"/dwi/*nii.gz"))

    niftis = glob.glob(outSesDir+"/*/*nii.gz")
    jsons = glob.glob(outSesDir+"/*/*json")

    # Get the files to remove, rename, or move to another directory
    to_remove, to_change = verify_output(subj, sess, niftis, jsons)
    master = (to_remove + to_change)
  
  
    # Execute the changes based on the list created above
    made_changes = correction(subj, sess, master, probs)

    subj = f"sub-{subj}"
    sess = f"ses-{sess}"
    

    # Write heudiconv's performance review to dataframe for initial run
    if relaunch is False:
        if actual_series_n < target_n:
            update_dict = create_dict(subj, sess, target_n, actual_series_n, actual_anat_n, anat_n, actual_dwi_n, dwi_n, 'missed')
            resubmit = True
        else:
            update_dict = create_dict(subj, sess, target_n, actual_series_n, actual_anat_n, anat_n, actual_dwi_n, dwi_n, 'success')
        
        outDf.loc[len(outDf)] = update_dict
              
    # Write heudiconv's performance review to dataframe for relaunched run
    else:
        
        if performance == 'fail' or actual_series_n > target_n:
            outDf.loc[(outDf.subject_id == subj)&(outDf.session_id == sess),'final_run'] = 'incomplete'
            
        else:
            outDf.loc[(outDf.subject_id == subj)&(outDf.session_id == sess),'final_run'] = 'success'
            resubmit = True

    return resubmit, made_changes


def create_dict(subj, sess, target_n, actual_series_n, actual_anat_n, anat_n, actual_dwi_n, dwi_n, performance):
    dct = {'subject_id' : subj,
        'session_id' : sess,
        'expected_output' : target_n,
        'actual_output' : actual_series_n,
        'expected_anat_scans' : anat_n,
        'actual_anat_scans': actual_anat_n,
        'expected_dwi_scans' : dwi_n,
        'actual_dwi_scans': actual_dwi_n,
        'first_run' : performance}
    return dct


def searchTxtFile(files, subj, sess):

    '''
    Reads in a text file to see if it contains a subject and session name

    Args:
        files: File to read
        subj: Subject name
        sess: Session name
    
    Return
        file: File containing the subject and session name
    '''
    
    found = []

    for file in files:
        with open(file, 'r') as f:
            for line in file:
                if subj in line and sess in line:
                    found.append(file)

        if len(found) == 1:
            return file

       
def fetch_answers(file, expected_series_n, relaunch=False):

    '''
    Read in slurm log of heudiconv's processing to evaluate it's performance

    Args:
        file: List of text files containing output of heudiconv processing
        relaunch: True/False if initial run of heudiconv or relaunch

    Return:
        actual_series_n: Number of expected output niftis based on dicoms that passed inclusion criteria
        anat_n: Number of expected anatomical niftis
        performance: success/fail for heudiconv's overall processing status
    '''

    actual_series_n = 0
    expected_scans_n = 0
    anat_n = 0
    dwi_n = 0
    performance = None

    print(file)
    try:
        file = file[0]
    except:
        print("Could not locate slurm logs with heudiconv processing.")
        sys.exit(1)
        
    

    with open(file, 'r') as f:
        for line in f:
            if relaunch is False:
                pattern = r'Generated sequence info for [0-9]+ studies with [0-9]+ entries total'
            else:
                pattern = r'Processing \b(?![01]\b)\d{1,4}\b pre-sorted seqinfo entries|Generated sequence info for 1 studies with 1 entries total'
            
            actual_series_n_line = re.findall(pattern, line) 
            performance_line = re.findall(r'INFO: PROCESSING DONE|json.decoder.JSONDecodeError|IndexError: list index out of range|RuntimeError: was asked to move', line)
            anat_n_line = re.findall(r'{bids_subject_session_dir}\/anat', line)
            expected_scans_n_line = re.findall(r'{bids_subject_session_dir}\/(anat|dwi)', line)
            dwi_n_line = re.findall(r'{bids_subject_session_dir}\/dwi', line)

            # Count the number of scans that pass through the heuristic that should be expected in the BIDS directory
            if expected_scans_n_line:
                expected_scans_n = expected_scans_n + 1
                
                # Count the number of anatomical scans that pass through the heuristic that should be expected in the BIDS directory
                if 'anat' in expected_scans_n_line:
                    anat_n = anat_n + 1

                # Count the number of anatomical scans that pass through the heuristic that should be expected in the BIDS directory
                if dwi_n_line: 
                    dwi_n = dwi_n + 1

            # Get the number of series that heudiconv registered
            if actual_series_n_line and relaunch is False:
                actual_series_n = int(actual_series_n_line[0].split("with")[1].split(" entries")[0])
                
            # Get the number of series that heudiconv registered after relaunching
            if actual_series_n_line and relaunch is True:
                if "pre-sorted" in actual_series_n_line[0]:
                    actual_series_n = actual_series_n + int(actual_series_n_line[0].split("Processing ")[1].split(" pre-sorted")[0])
                else:
                    actual_series_n = actual_series_n + int(actual_series_n_line[0].split("with ")[1].split(" entries")[0])
            
            # Get the overall performance of heudiconv 
            if performance_line:
                performance = 'fail' if 'fail' in performance_line[0] or 'Error' in performance_line[0] else "success"

    print(f"Heudiconv performance: {performance},\n number of series heudiconv registered: {expected_series_n},\n number of series heudiconv registered: {actual_series_n}.")
    if expected_series_n != actual_series_n:
        print("Relaunch needed due to missing series from heudiconv run")
    else:
        pass
    print(f"Number of expected niftis in BIDS: {expected_scans_n},\n\tanatomical: {anat_n}\n\tdiffusion: {dwi_n}")
    return actual_series_n, expected_scans_n, anat_n, dwi_n, performance


def iterate_run(inFn, lst, changes=False):
    
    '''
    Increase the run number of a file to prevent overwriting files with same acquisition parameters.

    Args:
        niftiInFn: Input filename without suffix (i.e. no nii.gz, json, etc.)
        lst: List of filenames
        changes: List of changes to make
    '''
    print("\t iterating run . . .")
    toMatch = re.sub(r'w[0-9]+','w', op.basename(inFn))
    if 'heudiconv' in toMatch:
        toMatch = toMatch.split("_heudiconv")[0] + re.findall(r'.nii.gz|.json|.bvec|.bval', toMatch)[0]
    
    matches = [i for i in lst if(i == toMatch)]
    iteration = int(len(matches)) + 1

    runNum = int(re.findall(r'run-\d+',inFn)[0].split("run-")[1].lstrip("0")) + iteration
    run = "run-"+ str(runNum).zfill(3)
    outFn = re.sub(r'run-\d+', run, inFn)
    lst.append(outFn)
    
    if changes == False:

        return outFn
        
    else:
        
        changes.append([inFn, outFn])

        for suffix in ['json','bvec','bval']:
            inAccFile = inFn.replace('nii.gz', suffix)
            if op.exists(inAccFile):
                outAccFile = outFn.replace("nii.gz", suffix)
                changes.append([inAccFile, outAccFile])
        

def correct_fn(map, subj, sess, outSesDir):

    '''
    Temporary function to combat CuBIDS only retaining some BIDS keys.
    Adds back in keys CuBIDS removed (sample, ce, mt, echo, part)

    Args:
        map: Bash script containing commands to change filenames
        subj: Subject name
        sess: Session name
        outSesDir: output subject session directory
    '''

    oldFn = None
    newFn = None
    substr = ["sample-fetal","ce-gad","mt-on","echo","part"]
    filenames = []
    changes = []
    pattern = f'sub-{subj}_ses-{sess}'
    made_changes = False
    
    niftis = glob.glob(outSesDir+'/*/*nii.gz')
    jsons = glob.glob(outSesDir+'/*/*json')

    with open(map, 'r') as file:
        for line in file:
            line = line.replace("\n","")
            if "git" in line and pattern in line:
                git, cmd, oldFn, newFn = line.split(" ")
            elif "mv" in line and pattern in line:
                cmd, oldFn, newFn = line.split(" ")
    

            if newFn in niftis:
                
                basename = newFn.rsplit("/",1)[1].split(".nii")[0]
                pairedJson = [i for i in jsons if (basename in i)][0]

                res = [word for word in substr if word in oldFn]
                
                if res:
                    name_parts = [
                        None if not "sample" in oldFn else "sample-fetal",
                        re.findall(r'acq-[A-z0-9]+_', newFn)[0].split("_")[0],
                        None if not "ce" in oldFn else "ce-gad",
                        None if not "rec" in oldFn else re.findall(r'rec-[A-Z0-9]+', oldFn)[0],
                        None if not "inv" in oldFn else re.findall(r'inv-[0-9]', oldFn)[0],
                        re.findall(r'run-[0-9]+', newFn)[0],
                        None if not "echo" in oldFn else re.findall(r'echo-[0-9]', oldFn)[0],
                        None if not "flip" in oldFn else re.findall(r'flip-[0-9]', oldFn)[0],
                        None if not "mt-on" in oldFn else "mt-on",
                        None if not "part" in oldFn else re.findall(r'part-[a-z]+', oldFn)[0],
                        newFn.rsplit("_",1)[1].split('.nii.gz')[0],                             # T1w, T2w, FLAIR, etc...
                    ]
                    
                    new_basename = pattern + "_" + "_".join(filter(bool, name_parts)) 
                    corrected_n = newFn.replace(basename, new_basename)
                    
                    '''
                    Check if the new name for the file already exists to prevent overwriting
                    If it exists, increase the run number
                    '''
                    if corrected_n in filenames:
                        iterate_run(corrected_n, filenames, changes)
                    else:
                        filenames.append(corrected_n)
                        changes.append([newFn, corrected_n])
                        corrected_j = corrected_n.replace("nii.gz","json")
                        changes.append([pairedJson, corrected_j])
                
                else:
                    if newFn not in filenames:
                        filenames.append(newFn)
                    else:
                        '''
                        Check if the new name for the file already exists to prevent overwriting
                        If it exists, increase the run number
                        '''
                        iterate_run(newFn, filenames, changes)
        
        if len(changes) > 0:
            made_changes = correction(subj=subj, sess=sess, lst=changes)
    return made_changes


def update_manifest(outSesDir):

    '''
    Update BIDS session manifest

    Args:
        outSesDir: output subject session directory
    '''

    manifestFn = glob.glob(outSesDir+"/*scans.tsv")[0]
    niftis = [i.split("BIDS/")[1] for i in glob.glob(outSesDir+"/*/*nii.gz")]
    df = pd.DataFrame(data = niftis, columns = ['filename'])
    df['acq_time'] = 'n/a'
    df.to_csv(manifestFn, index = False, sep = '\t')


def clean_up(fullpath, base, subj=None, sess=None):
    
    '''
    Remove heudiconv processing files

    Args:
        subj: Subject id
        sess: Session id
        base: Full path to the subjects' session BIDS subdirectory
        fullpath: Full path to the heudiconv hidden directory for a specific subject and session
    '''

    if op.exists(base):
      if subj and sess:
        scans = op.join(base, '*scans.tsv')
        if op.isfile(scans):
          os.system(f"rm {scans}")
        os.system(f"rm -rf {base}/*/*nii.gz")
        os.system(f"rm -rf {base}/*/*json")

    if op.exists(fullpath):
        dirsToRemove = [i for i in glob.glob(fullpath+"/*") if (op.isdir(i) and '.git' not in i)]
        for dirR in dirsToRemove:
            if op.exists(dirR):
                os.system(f"rm -rf {dirR}")


def copy_dicoms(sess, dcmDir, inSesDir): #copy_dicoms(sess, dcmDir, dcmName, inSesDir):
    
    tracker = []
    #dcmDir = dcmDir.rstrip("/",1)[0]
    if not op.exists(inSesDir):
        os.mkdir(inSesDir)

    # Copy session dicoms from source to input directory
    print()
    print(f"COPYING DICOMS FOR {sess} FROM {dcmDir}")
    util = f"cp -R {dcmDir}/* {inSesDir}"
    os.system(util)
    count = 1

    for content in glob.glob(f"{inSesDir}/*"):
        if '.zip' in content:
            
            extract_dir = content.replace(".zip","")
            print()
            print(f"EXTRACTING {content} to {extract_dir}")
            shutil.unpack_archive(content, extract_dir)

            dicom_set_dirs = glob.glob(f"{extract_dir}/**/*dcm",recursive=True)
            for dicom_set in dicom_set_dirs:
                dicom_set_dir = dicom_set.rsplit("/",1)[0]
                dicom_set_dirName = op.basename(dicom_set_dir)
                if dicom_set_dirName not in tracker:
                    tracker.append(dicom_set_dirName)
                    new_name = dicom_set_dir.replace(dicom_set_dirName,f"set_{count}")
                    print(f"RENAMING {dicom_set_dir} to {new_name}")
                    os.rename(dicom_set_dir, new_name)
                    count = count+1
                else:
                    continue
    

def run_heudiconv(fstr, outDir, subj, sess, dcmconfig, heuristic, logOutputDir, overwrite=False):

    print()
    
    if op.exists(op.join(outDir,'heudiconv.lock')):
        try:
            os.remove(op.join(outDir, 'heudiconv.lock'))
        except:
            pass

    heuristic = '/app/scripts/heuristic.py'

    if overwrite is False:
        
        print(f"RUNNING HEUDICONV FOR {sess}") 
        logFn = op.join(logOutputDir, f"bidsify_{subj}_{sess}.txt")
        util = f"singularity run -e -C -B {bgd} --env 'HEUDICONV_FILELOCK_TIMEOUT=45' {simg} heudiconv --files {fstr} -o {outDir} -f {heuristic} -s {subj} -ss {sess} -c dcm2niix --dcmconfig {dcmconfig} -b 2>&1 | tee {logFn}"
    else:
        print(f"REDOING HEUDICONV FOR {sess}") 
        logFn = op.join(logOutputDir, f"bidsify_redo_{subj}_{sess}.txt")
        util = f"singularity run -e -C -B {bgd} --env 'HEUDICONV_FILELOCK_TIMEOUT=45' {simg} heudiconv --files {fstr} -o {outDir} -f {heuristic} -s {subj} -ss {sess} -c dcm2niix --dcmconfig {dcmconfig} -b notop --overwrite 2>&1 | tee -a {logFn}"
    print(util)
    os.system(util)


def remove_false_scans(dir, wrkDir):

    print("CHECKING FOR FALSE SCANS")
    util = f"singularity run -e -C -B {dir} -B {wrkDir} {simg} python /app/scripts/update-files.py {dir} nifti"
    os.system(util)


def remove_run_tracker(subj, sess):
    # Remove the run tracker
    proc = sess.split("procId")[0]
    util = "find /mnt/isilon/bgdlab_processing/Data/tmp_heudiconv -name" + f" '{subj}_{proc}_*txt'" +" -exec rm {} \;"
    print("Removing run tracker")
    print(util)
    os.system("find /mnt/isilon/bgdlab_processing/Data/tmp_heudiconv -name" + f" '{subj}_{proc}_*txt'" +" -exec rm {} \;")


def main():

  '''
  Main function of the script
  '''

  parser = argparse.ArgumentParser()
  parser.add_argument("work_dir")
  parser.add_argument("-d", "--dicom-paths", 
                      help = "Full path to the directory(ies) containing the session dicoms.",
                      nargs = '+', 
                      required = False)
  '''parser.add_argument("-d", "--dicom-basepath", 
                      help = "Full path to the directory containing the session dicoms.",
                      required = False)
  parser.add_argument("-dd", "--dicom-dirnames", 
                      nargs = '+', 
                      help = "Directory names of the session dicoms.",
                      required = False)'''
  parser.add_argument("-i", "--input-directory", 
                      help = "Full path to the subject level directory to hold the dicom files",
                      required = False)
  parser.add_argument("-o", "--output-directory",
                        help = "Full path to the output directory for the bidsified data",
                        required = True)
  parser.add_argument("-s", "--subject", 
                      help="Subject ID",
                      required=True)
  parser.add_argument("-ss", "--sessions", 
                      nargs = '+', 
                      help = "Session ID or list of sessions ids",
                      required=True)
  parser.add_argument("-pc", "--post-cubids", 
                      action = 'store_true', 
                      help="Post-CuBIDS udpates")
  parser.add_argument("-mm", "--map", 
                      help = "CuBIDS bash script to change filenames",
                      required = False)

  args = parser.parse_args()
  wrkDir = args.work_dir
  #dcmDir = args.dicom_basepath
  #dcmDirNames = args.dicom_dirnames
  dcmPaths = args.dicom_paths
  srcDir = args.input_directory
  outDir = args.output_directory
  cubids = args.post_cubids
  subj = args.subject.strip()
  sessions = args.sessions
  

  # Pre CuBIDS organization: --------------------------------------------
  if cubids is False:
    outDf = pd.DataFrame(data = [], 
                        columns = ['subject_id','session_id',
                                    'expected_output',"actual_output",
                                    'expected_anat_scans','actual_anat_scans',
                                    'expected_dwi_scans','actual_dwi_scans',
                                    'first_run'])
    probs = pd.DataFrame(data = [],
                        columns = ['subject_id','session_id','path'])
    
    n_sess = len(sessions)
    heuristic = op.join(wrkDir, 'heuristic.py')
    logOutputDir = op.join(wrkDir, 'logs', 'heudiconv')
    
    print(f"\nCURATING SUBJECT {subj} DATA WITH {n_sess} SESSIONS")
    ################## COPY DICOMS ##################
    inSubjDir = op.join(srcDir, subj)
    outSubjDir = op.join(outDir, f"sub-{subj}")

    

    src = dl.Dataset(srcDir) 

    srcSub = dl.Dataset(inSubjDir) 
    if not srcSub.is_installed():

        print("CREATING SOURCE SUBJECT DATASET")  
        srcSub.create(cfg_proc='text2git')
    
    
    ################## PREP ##################
    for sess, dcmPath in zip(sessions, dcmPaths):
        
        print("----------------------------------------------------------------")
        inSesDir = op.join(inSubjDir, sess)
        srcSes = dl.Dataset(inSesDir)
        if not srcSes.is_installed():
            srcSes.create(cfg_proc='text2git')
        
        # Copy session dicoms from source to input directory
        copy_dicoms(sess, dcmPath, inSesDir)
        print("\nSAVING DICOM SESSION DATASET")
        srcSes.save(path=inSesDir, message="Added dicoms", jobs='auto', result_renderer = 'disabled')
        
        # Get the dicoms
        dcmFns = glob.glob(inSesDir+"/**/*dcm", recursive=True)

        srcSes.unlock(result_renderer = 'disabled')

        print("\nCHECKING IF SESSION HAS DICOMS IN EXPECTED FORMAT")
        if '.dcm' in dcmFns[0]:
            pass
        else:
            print(f"This sessions appears to not have *dcm files, check {inSesDir} for filetype.")
            sys.exit(1)
        
        # Filter and decompress dicoms in the container
        print(f"\nFILTERING DICOMS")
        
        util = f'singularity run -e -C -B {inSesDir} -B {wrkDir} {simg} python /app/scripts/update-files.py {inSesDir} dicom'
        result = subprocess.run(util.split(), capture_output=True, text=True)
        if result.returncode == 0:
            print(result.stdout)
        else:
            print_cmd(result)
            print(result.stderr)
            sys.exit(1)

        print("\nSAVING DICOM SESSION DATASET")
        srcSes.save(path=inSesDir, message='Edited metadata and decompressed dicoms.', jobs='auto', result_renderer = 'disabled')

    print("\nSAVING DICOM SUBJECT DATASET")
    srcSub.save(path=inSubjDir, message='Edited dicoms for sessions', jobs='auto', result_renderer = 'disabled')
    
    bidsSub = dl.Dataset(outSubjDir)
    if not bidsSub.is_installed():
        bidsSub.create(cfg_proc='text2git')

    ################## HEUDICONV RUNS ##################
    for sess in sessions:
        
        print("----------------------------------------------------------------")
        inSesDir = op.join(inSubjDir, sess)
        srcSes = dl.Dataset(inSesDir)
        outSesDir = op.join(outSubjDir, f"ses-{sess}")
        bidsSes = dl.Dataset(outSesDir) 
        if not bidsSes.is_installed():
            bidsSes.create(cfg_proc='text2git')

        srcSes.unlock(result_renderer = 'disabled')

        # Count subdirectories
        dirs = os.listdir(inSesDir)
        
        oneDcm = glob.glob(f"{inSesDir}/**/*dcm", recursive=True)[0]
        if 'set' not in oneDcm:
            dirsWithDcms = [i for i in dirs if ('GCP' in i or '1.' in i)]
            if len(dirsWithDcms) == 1:
                try:
                    files = glob.glob(f"{inSesDir}/GCP*/*")[0]
                except:
                    files = glob.glob(f"{inSesDir}/1.*")[0]
        else:
            # Find the parent directory with the dicom set subdirectories
            try:
                files = oneDcm.rsplit("/",2)[0]
            except:
                sys.exit(1)
        
        
        # Create the custom dcmconfig.json file
        dcmconfig = create_config(inSesDir)

        # Run heudiconv
        run_heudiconv(files, outDir, subj, sess, dcmconfig, heuristic, logOutputDir)
        
        print("\nSAVING BIDS SESSION DATASET")
        bidsSes.save(path=outSesDir, message='Initial heudionv run', jobs='auto', result_renderer = 'disabled')

        # Check output
        bidsSes.unlock(result_renderer = 'disabled')

        print(f"\nCHECKING OUTPUT ses-{sess}")
        remove_false_scans(outSesDir, wrkDir)
        redo_heudiconv, made_changes = check_output(wrkDir, subj, sess, files, outSesDir, probs, outDf)
        
        
        if made_changes is True or is_clean(outSubjDir) is False: 
            print("\nSAVING BIDS SESSION DATASET")
            bidsSes.save(path=outSesDir, message="Removing unpaired niftis/jsons and non-brain scans. Moving errant DWI scans and confusing anat scans to problematic_scans dir. Name correction for multi-echo scans.",
                            jobs='auto',
                            result_renderer = 'disabled')

        ''' 
        Redo heudiconv series wise to get around heudiconv stopping 
        once it encounters an error when running session wise.
        '''  
        if redo_heudiconv is True:
            
            srcSes.unlock(result_renderer = 'disabled')
            bidsSes.unlock(result_renderer = 'disabled')

            remove_run_tracker(subj, sess)

            # Remove the session scan manifest, modality directories, 
            # and anything in the problematics scans directory
            clean_up(outSesDir, outSesDir, subj, sess)
            
            dcmDirs = [i.rsplit("/",1)[0] for i in glob.glob(f"{inSesDir}/**/*.dcm", recursive=True)]
            dcmDirs = list(set(dcmDirs))
            
            for series in dcmDirs:        

                # Remove heudiconv processing tables and heudiconv lock if it exists
                clean_up(op.join(outDir, '.heudiconv', subj), outSesDir)
                run_heudiconv(series, outDir, subj, sess, dcmconfig, heuristic, logOutputDir, True)

            print("\nSAVING BIDS SESSION DATASET")
            bidsSes.save(path=outSesDir,message="Heudiconv redo", jobs='auto', result_renderer = 'disabled')
            bidsSes.unlock(result_renderer = 'disabled')

            # Check output
            print(f"\nCHECKING OUTPUT ses-{sess}")
            remove_false_scans(outSesDir, wrkDir)
            redo_heudiconv, made_changes = check_output(wrkDir, subj, sess, files, outSesDir, probs, outDf, True)
            if made_changes is True: 
                bidsSes.save(path=outSesDir, message="After heudiconv relaunch: removing unpaired niftis/jsons and non-brain scans. Moving errant DWI scans and confusing anat scans to problematic_scans dir. Name correction for multi-echo scans.",
                                jobs='auto',
                                result_renderer = 'disabled')

            print(f"\nCHANGING FILE PERMISSIONS")
            os.system(f"chmod +rw -R {outSesDir}")

            print("\nSAVING SESSION SUBJECT DATASET")
            bidsSes.save(path=outSesDir, message="Changing file permissions", jobs='auto', result_renderer = 'disabled')

        remove_run_tracker(subj, sess)

        print("\nSAVING DICOM SESSION DATASET")
        # print(is_clean(inSesDir))
        srcSes.save(path=inSesDir, message='Added dcmconfig', jobs='auto', result_renderer = 'disabled')

    print("\nSAVING DICOM SUBJECT DATASET")
    print(is_clean(inSubjDir))
    srcSub.save(path=inSubjDir, jobs='auto', result_renderer = 'disabled')

    
    
                    
    print("\nSAVING BIDS SUBJECT DATASET")
    print(is_clean(outSubjDir))
    bidsSub.save(path=outSubjDir, message='Session curation', jobs='auto', result_renderer = 'disabled')
          
          
      
    ################## WRITE OUTPUT ##################
    print()
    print(f"WRITING OUTPUT")

    problematicScansDir = op.join(outDir.rsplit("/",1)[0],'problematic_scans')
    if not probs.empty and op.exists(problematicScansDir):  
        probs_path = op.join(problematicScansDir,f"{subj}_{sess}_problematic_scans.tsv")
        probs.drop_duplicates(subset=["path"]).to_csv(probs_path,sep="\t",index=False)

    outDfFn = op.join(wrkDir,"logs","heudiconv",f"{subj}_heudiconv_perform.tsv")
    outDf.to_csv(outDfFn, index=False, sep = "\t")

      


  # Post CuBIDS organization: --------------------------------------------
  else:

    fnMap = args.map
    outSubjDir = op.join(outDir, f"sub-{subj}")
      

    bidsSub = dl.Dataset(outSubjDir)
    
    ################## MAKING UPDATES ##################
    for sess in sessions:
        outSesDir = op.join(outSubjDir, f"ses-{sess}")
        bidsSes = dl.Dataset(outSesDir)

        if not op.exists(outSesDir):
            print(f"Directory {outSesDir} does not exist. Either because a problem with heudiconv or because no data passed heuristic filter")
            sys.exit(1)
        else:
            pass
        
        bidsSes.unlock(result_renderer = 'disabled')
        
        print("\nUPDATING SESSIONS' SCAN MANIFEST")
        sess_manifest = op.join(outSesDir, f'sub-{subj}_ses-{sess}_scans.tsv')
        bidsSes.unlock(path=sess_manifest)
        update_manifest(outSesDir)
        bidsSes.save(path = outSesDir,
                        message = 'Updating session scan manifests.', 
                        jobs='auto',
                        result_renderer = 'disabled')
        
    bidsSub.save(path = outSubjDir,
                    message = "Post CuBIDS updates.", 
                    jobs='auto',
                    result_renderer = 'disabled')

    
if __name__ == "__main__":
    main()
