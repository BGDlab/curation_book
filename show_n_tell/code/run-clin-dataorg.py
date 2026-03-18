import argparse
import sys
import os
import os.path as op
import pandas as pd
import re
import glob
import subprocess
from datetime import datetime
import datalad.api as dl

""" 

Usage

    If subject ID pattern isn't HM[0-9A-Z]+, use -p <regex> 

    Setup and organize at subject level 

        python run-clin-dataorg.py <path-to-working-directory> bidsify -b <path-to-base-directory> -d <path-to-input-dicoms> -n <name-of-dataset>

    Run CuBIDS (prep and apply)

        python run-clin-dataorg.py <path-to-working-directory> cubids -b <path-to-base-directory> -n <name-of-dataset> 

    Post-CuBIDS corrections

        python run-clin-dataorg.py \
            <path-to-working-directory> \
            post-cubids \
            -b <path-to-base-directory> \
            -n <name-of-dataset> \
            -d <path-to-input-dicoms> \
            -m <path-to-delivery_manifest.csv> \
            -mm <path-to-cubids-run-script.sh> 

    Check job status from step 1:

        python run-clin-dataorg.py <path-to-working-directory> bidsify -j

    Check job status from step 3:

        python run-clin-dataorg.py <path-to-working-directory> post-cubids -j 
    
    Combine batched CuBIDS apply output:
    
        python -u run-clin-dataorg.py <path-to-working-directory> combine -b <path-to-CuBIDS-batch>

"""

# Some global things
today = datetime.today().strftime('%Y-%m-%d')
container_url = 'docker://dabrielz/clin-dataorg:0.0.2'
simg = '/mnt/isilon/bgdlab_processing/code/containers/clin-dataorg-0.0.2.sif'
bgd = '/mnt/isilon/bgdlab_processing/'
fs22q11 = '/mnt/isilon/22q11/projects/'

    
def cli():
    ''' 
    CLI for container 
    '''
    parser = argparse.ArgumentParser(description='Clinical data organization')

    parser.add_argument('work_dir',
        help='Path to the workding directory where the repo is stored.')

    parser.add_argument('cmd',
        help='Command to run',
        choices=['setup', 'bidsify', 'cubids-prep',
                 'cubids', 'cubids-apply', 'post-cubids',
                 'annotate', 'combine'])

    parser.add_argument('-b', '--base',
        help='Path to the base directory, default is /mnt/isilon/bgdlab_processing/Data/SLIP.',
        default='/mnt/isilon/bgdlab_processing/Data/SLIP',
        required=False)

    parser.add_argument('-n' ,'--dataset_name',
        help='The name of the dataset, e.g. slip_2023_03 or 2023-07_genetics_patient_imaging.',
        default='',
        required=False)

    parser.add_argument('-d', '--dicoms',
        help='Path to the delivered dicoms, in the 22q11 fileshare.',
        nargs = '+', 
        default='',
        required=False)
    
    parser.add_argument('-m','--manifest',
        help='Path to the request manifest csv file.',
        default='',
        required=False)
    
    parser.add_argument('-p', '--pattern',
        help='The regex pattern of the subject label, default is HM[A-Z0-9]+.',
        default='HM[A-Z0-9]+',
        required=False)
    
    parser.add_argument('-mm','--cubids-full-cmd',
        help='Path to the bash script created by CuBIDS to change filenames.',
        required=False)

    parser.add_argument('-j', '--job-summary',
        help='Print job summary for indicated stage.',
        action='store_true')
    
    parser.add_argument('-test',
        help='Run code tests',
        action='store_true')
    
    parser.add_argument('-r', '--relaunch',
        help='Relaunch failed jobs',
        action='store_true')
    
    parser.add_argument('--partial',
        help='Limit number of concurrent running array jobs',
        action='store_true')

    return parser


def print_cmd(result):

    '''
    Print the command used for a subprocess run

    Args:
        result: Output of subprocess run
    '''
    fstr =" ".join(result.args).replace("'","")
    print(fstr)


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


def data_setup(wrkDir, manifest, batch, pattern, outFn):
    '''
    Prep data for curation: 
    1) check for duplicate delivery (e.g. proc ids were included in a previous delivery)
    2) Verify no multi-uuid

    Args:
        manifest: Delivery manifest table
        batch: Full path to the directory(ies) containing the raw dicoms 
        pattern: Subject ID pattern to look for
        outFn: Full path to output filename
    
    Return:
        subprocess output from running data-prep.py
    '''

    print()
    print()
    print("===== PREPARING DATA =====")
    print()
    
    util = f"singularity run -e -C -B {bgd} -B {fs22q11} {simg} python /app/scripts/data-prep.py {wrkDir} -i {batch} -t {manifest} -s {pattern} -o {outFn}"
    result = subprocess.run(util.split(), capture_output=True,  text=True)

    
    return result


def data_curate(wrkDir, config, bidsDir, partial):
    
    '''
    Launch an array of jobs to:
    1) copy dicoms from the input source
    2) remove problematic metadata tags and decompress dicoms
    3) run heudiconv 
    4) check heudiconv output and remove problematic scans  
    at the subject level

    Args:
        config: Table of subjects with information used to launch 
            jobs to curate data
        srcDir: Full path to the dicom store (sourcedata directory)
        bidsDir: Full path to the output directory (BIDS directory)
        partial: True/False to limit number of concurrent running jobs
    
    Return:
        result: subprocess run output
        n: Number of jobs launched
        jobid: Array job id

    '''

    print()
    print()
    print("===== CURATING DATA =====")
    print()
    print(config)
    configDf = pd.read_csv(config)

    n = configDf.pat_id.nunique()
    n_jobs = n + 1
    if partial:
        strArray = f"2-{n_jobs}%200"
    else:
        strArray = f"2-{n_jobs}"
    

    curLogOutput = op.join(wrkDir,"logs","curation")
    heudLogOutput = op.join(wrkDir, 'logs', 'heudiconv')
    if not op.exists(curLogOutput):
        os.makedirs(curLogOutput)
    if not op.exists(heudLogOutput):
        os.mkdir(heudLogOutput)


    arrayScript = op.join(wrkDir, "array-launcher.sh")
    dataset = bidsDir.replace("/BIDS","")
    util = f"sbatch -a {strArray} --mem-per-cpu 5G -o {curLogOutput}/curate-%A_%a.out --time 24:00:00 {arrayScript} {wrkDir} {dataset} False {config}"
    result = subprocess.run(util.split(), capture_output=True, text=True)
    print_cmd(result)
    
    if result.returncode == 0:    
        jobid = str(result.stdout.split("Submitted batch job ")[1]).split("\n")[0]
        print(f"LAUNCHED {n} JOBS. Job: {jobid}")     
        jobs = [f"{jobid}_{i}" for i in range(2, (2+n))] 
        jobsummary = pd.DataFrame({'pat_id':configDf.pat_id.unique(),
                                   'jobid':jobs})
        jobSumOut = op.join(wrkDir, 'curation_job_summary.tsv')
        if not op.exists(jobSumOut):
            jobsummary.to_csv(jobSumOut, index = False, sep = '\t')
        else:
            prevDf = pd.read_csv(jobSumOut, sep = '\t')
            prevDf['origin'] = 'ORIGINAL'
            jobsummary['origin'] = 'RELAUNCH'
            merged = pd.concat([prevDf, jobsummary])
            merged.to_csv(jobSumOut, index = False, sep = '\t')


        
    else:
        print('FAILED TO LAUNCH JOBS FOR SUBJECT LEVEL CURATION')
        print(result.stderr)
    
    return result, n, jobid


def get_job_summary(path):

    '''
    Output a summary of the past jobs submitted.

    Args:
        path: Full path to the table containing job ids to check the status of 
    '''

    df = pd.DataFrame(pd.read_csv(path, sep = '\t'))
    
    for idx, row in df.iterrows():
        jobid = row['jobid']
        result = subprocess.run(['sacct', f"--jobs={jobid}", "--format=jobname,state,submit,end","-p"], capture_output = True, text = True)
        jobstatus = re.findall('COMPLETED|FAILED|RUNNING|PENDING|OUT_OF_MEMORY|TIMEOUT|CANCELLED', result.stdout)
        if not jobstatus:
            jobstatus = 'NOT LAUNCHED'
            df.loc[idx, 'status'] = jobstatus
        else:
            jobstatus = jobstatus[0]
            df.loc[idx, 'status'] = jobstatus
            try:
                start = re.findall(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}',result.stdout)[0]
                df.loc[idx, 'start'] = pd.to_datetime(start)
            except:
                pass
            try:
                end = re.findall(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}',result.stdout)[1]
                df.loc[idx, 'end'] = pd.to_datetime(end)
            except:
                pass

    if 'end' in df.columns and 'start' in df.columns:
        df_fil = df[(df['end'].notna())&(df['start'].notna())]
        diff = df_fil.end.max()-df_fil.start.min()
        print(f"Array time elapsed: {diff}")
    df.to_csv(path, sep = "\t", index = False)
    
    
    if 'origin' not in df.columns:
        n = df.pat_id.nunique()
        success = df[df['status']=='COMPLETED'].pat_id.nunique()
        fail = df[df['status']=='FAILED'].pat_id.nunique()
        running = df[df['status']=='RUNNING'].pat_id.nunique()
        pending = df[df['status']=='PENDING'].pat_id.nunique()
        oom = df[df['status']=='OUT_OF_MEMORY'].pat_id.nunique()
        tm = df[df['status']=='TIMEOUT'].pat_id.nunique()
        canc = df[df['status']=='CANCELLED'].pat_id.nunique()
        na = df[df['status']=='NOT LAUNCHED'].pat_id.nunique()
        print(f"\nOut of {n} jobs:")
        message = [
            None if not na else "%d not yet launched" % na,
            None if not pending else "%d are still pending" % pending if pending > 1 else "1 is still pending",
            None if not running else "%d are still running" % running if running > 1 else "1 is still running",
            None if not oom else "%d stopped, OOM" % oom,
            None if not tm else "%d stopped, timed out" % tm,
            None if not canc else "%d were cancelled" % canc if canc > 1 else "1 was cancelled",
            None if not fail else "%d failed" % fail,
            None if not success else "%d succeeded" % success,
        ] 
        message = "\n\t".join(filter(bool, message)) 
        print(f"\t{message}")
    else:
        og = df[df['origin']=='ORIGINAL']
        success = og[og['status']=='COMPLETED'].pat_id.nunique()
        fail = og[og['status']=='FAILED'].pat_id.nunique()
        running = og[og['status']=='RUNNING'].pat_id.nunique()
        pending = og[og['status']=='PENDING'].pat_id.nunique()
        oom = og[og['status']=='OUT_OF_MEMORY'].pat_id.nunique()
        tm = og[og['status']=='TIMEOUT'].pat_id.nunique()
        canc = og[og['status']=='CANCELLED'].pat_id.nunique()
        na = og[og['status']=='NOT LAUNCHED'].pat_id.nunique()
        n = og.pat_id.nunique()
        print(f"\nOut of {n} jobs from original launch:")
        message = [
            None if not na else "%d not yet launched" % na,
            None if not pending else "%d are still pending" % pending if pending > 1 else "1 is still pending",
            None if not running else "%d are still running" % running if running > 1 else "1 is still running",
            None if not oom else "%d stopped, OOM" % oom,
            None if not tm else "%d stopped, timed out" % tm,
            None if not canc else "%d were cancelled" % canc if canc > 1 else "1 was cancelled",
            None if not fail else "%d failed" % fail,
            None if not success else "%d succeeded" % success,
        ] 
        message = "\n\t".join(filter(bool, message)) 
        print(f"\t{message}")

        og = df[df['origin']=='RELAUNCH']
        success = og[og['status']=='COMPLETED'].pat_id.nunique()
        fail = og[og['status']=='FAILED'].pat_id.nunique()
        running = og[og['status']=='RUNNING'].pat_id.nunique()
        pending = og[og['status']=='PENDING'].pat_id.nunique()
        oom = og[og['status']=='OUT_OF_MEMORY'].pat_id.nunique()
        tm = og[og['status']=='TIMEOUT'].pat_id.nunique()
        canc = og[og['status']=='CANCELLED'].pat_id.nunique()
        na = og[og['status']=='NOT LAUNCHED'].pat_id.nunique()
        n = og.pat_id.nunique()
        print(f"\nOut of {n} jobs from relaunch:")
        message = [
            None if not na else "%d not yet launched" % na,
            None if not pending else "%d are still pending" % pending if pending > 1 else "1 is still pending",
            None if not running else "%d are still running" % running if running > 1 else "1 is still running",
            None if not oom else "%d stopped, OOM" % oom,
            None if not tm else "%d stopped, timed out" % tm,
            None if not canc else "%d were cancelled" % canc if canc > 1 else "1 was cancelled",
            None if not fail else "%d failed" % fail,
            None if not success else "%d succeeded" % success,
        ] 
        message = "\n\t".join(filter(bool, message)) 
        print(f"\t{message}")


def cubids_prep(bidsDir):

    '''
    Run CuBIDS steps add-nifiti-info, group, and apply

    Args:
        bidsDir: Full path to the output directory (BIDS directory)
    '''
    
    print()
    print()
    print("===== CuBIDS Prep =====")
 
    cubids_config = '/app/scripts/config.yml'
    script = '/app/scripts/identifyFilesFromCubids.py'
    cubidsDir = op.join(bidsDir, "code", "CuBIDS")
    parameters ='/app/scripts/parameter_ranges.tsv'
    
    print("UNLOCKING")
    ds = dl.Dataset(bidsDir)
    ds.unlock(recursive=True, result_renderer='disabled')

    print("ADDING NIFTI INFO")
    util = f"singularity run -e -C -B {bgd} {simg} cubids add-nifti-info {bidsDir}"
    addN = subprocess.run(util.split(), capture_output=True, text=True)
    print_cmd(addN)

    if addN.returncode == 0:
        
        ds.save(path=bidsDir, recursive=True, message='Add nifti info.', jobs = 'auto', result_renderer='disabled')

        print("GROUPING")
        ds.unlock(recursive=True, result_renderer='disabled')
        util = f"singularity run -e -C -B {bgd} {simg} cubids group --config {cubids_config} {bidsDir} run01"
        grouping = subprocess.run(util.split(), capture_output=True, text=True)
        print_cmd(grouping)
       
        ds.save(path=bidsDir, message="CuBiDS group run01", result_renderer='disabled')

        if grouping.returncode == 0:
            
            print("SCAN IDENTIFICATION")
            ds.unlock(path=op.join(cubidsDir, 'run01_summary.tsv'))
            util = f"singularity run -e -C -B {bgd} {simg} python {script} -s {cubidsDir}/run01_summary.tsv -g {parameters} -d /app/scripts"
            identification = subprocess.run(util.split(), capture_output=True, text=True)  
            print_cmd(identification)
            
            if identification.returncode == 0:
                print(identification.stdout)

                ds.save(path=bidsDir, message="Scan identification run01", result_renderer='disabled')
 
            else:
                print("SCAN IDENTIFICATION FAILED")
                print(identification.stderr)   
        else:
                print("GROUPING FAILED")
                print(grouping.stderr) 
    else:
        print("FAILED TO ADD NIFTI INFO")
        print(addN.stderr)    

    
def cubids_apply(wrkDir, deliveryPath, bidsDir, dependency=False):
    
    '''
    Run cubids-apply in the containerized environment. 
    Run normally if less than 4000 scans, run 10 batches if more than 4000 scans.

    Args:
        wrkDir: Working directory where commands are launched from
        deliveryPath: Directory containing the BIDS dataset
        bidsDir: Path to the BIDS dataset
        dependency: SLURM ID of cubids-prep job to use as dependency
    '''
    
    print("===== CuBIDS Apply =====")
    
    ds = dl.Dataset(bidsDir)

    applyScript = op.join(wrkDir, 'cubids-apply.sh')
    codeDir = f"{bidsDir}/code/CuBIDS"
    batchCombine = op.join(wrkDir, 'batch-cubids-combine.sh')

    df = pd.read_csv(f"{codeDir}/run01_files.tsv", sep = '\t')
    print("Number of files:", len(df))

    if len(df) <= 4000:
        util = f"sbatch --mem-per-cpu 8G -o {wrkDir}/cubids_apply_%j.out --time 8-12:00 {applyScript} {wrkDir} {deliveryPath} False"
    
    else:
        
        # Batching the CuBIDS apply run so that the process goes faster.  

        if not op.exists(op.join(codeDir, 'batch')):
            os.mkdir(op.join(codeDir, 'batch'))

        df2 = pd.read_csv(f"{codeDir}/run01_summary_edited.tsv", sep = '\t')
        df2_filtered = df2[(df2['RenameEntitySet'].notna())]
        toRemove = False
        entity_sets = df2_filtered[(df2_filtered['MergeInto']!=0)]['RenameEntitySet'].unique()
        n_sets_tot = df2_filtered[(df2_filtered['MergeInto']!=0)]['RenameEntitySet'].nunique()
        print(f"Number of RenameEntitySet: {n_sets_tot}")
        
        # First need to sort the edited summary file to get all the rows with data to be removed.
        if 0 in df2.MergeInto.unique():
            toRemove = True
            df2_remove = df2_filtered[df2_filtered['MergeInto']==0]
            n_remove_sets_tot = df2_remove.RenameEntitySet.nunique()
            print(f"Number of RenameEntitySets to remove: {n_remove_sets_tot}")
        
        splits = list(split(entity_sets, 10))
        split_config = pd.DataFrame(columns=['split_id','summary_path','files_path','name'])
        split_id = 1
        
        while split_id <= 10:

            # Get the sets to include in each split
            sets = splits[split_id-1]

            # Split the tables
            split_sets = df2_filtered[df2_filtered['RenameEntitySet'].isin(sets)]
            params = split_sets.KeyParamGroup.tolist()
            n_params = split_sets.KeyParamGroup.nunique()
            split_files = df[df['KeyParamGroup'].isin(params)]
            n_files = split_files.FilePath.nunique()

            # Setup output
            split_sets_path = op.join(codeDir, 'batch', f'batch{split_id}_summary_edited.tsv')
            split_files_path = op.join(codeDir, 'batch', f'batch{split_id}_files.tsv')

            # Write split output
            split_sets.to_csv(split_sets_path, sep = '\t', index=False)
            split_files.to_csv(split_files_path, sep = '\t', index=False)

            # Update main tables
            df2.loc[df2['KeyParamGroup'].isin(params), 'batch'] = split_id
            df.loc[df['KeyParamGroup'].isin(params), 'batch'] = split_id

            # Update config
            split_config.loc[len(split_config)] = {"summary_path": split_sets_path,
                                                   "files_path": split_files_path,
                                                   "name": f"run02_batch{split_id}"}
                                
            print(f"batch {split_id} contains {n_params} ParamGroups and {n_files} files to rename.")
            split_id = split_id + 1
        
        if toRemove == True:
            
            # Setup output
            remove_sets_path = op.join(codeDir, 'batch', f'batch{split_id}_summary_edited.tsv')
            remove_files_path = op.join(codeDir, 'batch', f'batch{split_id}_files.tsv')
            
            # Get the sets to remove
            params = df2_remove.KeyParamGroup.tolist()
            df_remove = df[df['KeyParamGroup'].isin(params)]
            n_params = df2_remove.KeyParamGroup.nunique()
            n_files = df_remove.FilePath.nunique()

            # Write split output
            df2_remove.to_csv(remove_sets_path,  sep = '\t', index=False)
            df_remove.to_csv(remove_files_path,  sep = '\t', index=False)

            # Update main tables
            df2.loc[df2['KeyParamGroup'].isin(params), 'batch'] = split_id
            df.loc[df['KeyParamGroup'].isin(params), 'batch'] = split_id
            
            # Update config
            print(f"batch {split_id} contains {n_params} ParamGroups and {n_files} files to rename.")
            split_config.loc[len(split_config)] = {"summary_path": remove_sets_path,
                                                   "files_path": remove_files_path,
                                                   "name": f"run02_batch{split_id}"}

        # Save tables
        split_config.to_csv(f"{wrkDir}/batch_cubids_apply_config.csv", index=False)
        df2.to_csv(f"{codeDir}/run01_summary_edited.tsv", sep = '\t', index=False)
        df.to_csv(f"{codeDir}/run01_files.tsv", sep = '\t', index=False)

        array = f"2-{split_id}"
        util = f"sbatch -J batch_apply_%a -a {array} --mem-per-cpu 12G -o {wrkDir}/cubids_apply_%A-%a.out --time 5-12:00 {applyScript} {wrkDir} {deliveryPath} {wrkDir}/batch_cubids_apply_config.csv"
        print(util)

    if dependency is not False:
        util = util.replace('sbatch', f"sbatch -d afterok:{dependency}")

    cApply = subprocess.run(util.split(), capture_output=True, text=True)
    print_cmd(cApply)

    if cApply.returncode != 0:
        print(cApply.stderr)
    else:
        jobid = str(cApply.stdout.split("Submitted batch job ")[1]).split("\n")[0]
        print(f"job: {jobid}")

        # Launch wrap up job to run after all batches run
        if len(df) > 4000:
            
            util = f"sbatch -d afterok:{jobid} -J combine_batches --mem-per-cpu 12G -o {wrkDir}/cubids_batch_combine_%j.out --time 12:00 {batchCombine} {wrkDir} {deliveryPath}"
            wrapUp = subprocess.run(util.split(), capture_output=True, text=True)
            
            print_cmd(wrapUp)
            if wrapUp.returncode != 0:
                print(wrapUp.stderr)
            else:
                jobid = str(wrapUp.stdout.split("Submitted batch job ")[1]).split("\n")[0]
                print(f"job: {jobid}")

    
def post_cubids(wrkDir, config, fnMap, bidsDir):
    '''
    Perform post-CuBIDS updates to include:
    1) Correcting filenames since CuBIDS doesn't currently retains all BIDS keys
    2) Update scan.tsv for each session
    at the subject level
    '''
    
    print()
    print()
    print(f"===== POST CuBIDS UPDATES =====")


    configDf = pd.DataFrame(pd.read_csv(config))
    
    n = configDf.pat_id.nunique()
    n_jobs = n + 1
    if n_jobs < 100:
        strArray = f"2-{n_jobs}"
    else:
        strArray = f"2-{n_jobs}%100"

    logOutput = op.join(wrkDir,"logs","curation")
    if not op.exists(logOutput):
        os.system(f'mkdir -p {logOutput}')


    arrayScript = op.join(wrkDir, "array-launcher.sh")
    dataset = bidsDir.replace("/BIDS","")
    util = f"sbatch -a {strArray} --mem-per-cpu 2G -o {logOutput}/postCubids-%A_%a.out --time 2:00:00 {arrayScript} {wrkDir} {dataset} True {config} {fnMap}"

    result = subprocess.run(util.split(), capture_output=True, text=True)
    print_cmd(result)
    
    if result.returncode == 0:  
        jobid = str(result.stdout.split("Submitted batch job ")[1]).split("\n")[0]
        print(f"LAUNCHED {n} JOBS. Job: {jobid}")   
        jobs = [f"{jobid}_{i}" for i in range(2, (2+n))] 
        jobsummary = pd.DataFrame({'pat_id':configDf.pat_id.unique(),'jobid':jobs})
        jobsummary.to_csv(op.join(wrkDir, 'postCubids_job_summary.tsv'), index = False, sep = '\t')

            
    else:
        print('FAILED TO LAUNCH JOBS FOR SUBJECT LEVEL CURATION')
        print(result.stderr)
    
    
    return result


def annotate(wrkDir, bidsDir, base, name, batch, manifest):
    
    """
    
    Annotate dataset, creating the README and accessory BIDS files.

    """
    
    util = f"singularity run -e -C -B {bgd} {simg} python -u /app/scripts/annotate.py {wrkDir} -b {base} -n {name} -d {batch} -m {manifest}"

    result = subprocess.run(util.split(), capture_output=True, text=True)
    print_cmd(result)
    if result.returncode == 0:
        print(result.stdout)
        
        
        ds = dl.Dataset(bidsDir)
        ds.save(path=bidsDir, message = "Adding participants files and README")

    else:
        print(result.stderr)
        sys.exit(1)


def concat_tables(tbls, output, bash=False): 
    dfList = []
    fstr = ["cat "]

    if bash == True:
        for tbl in tbls:
            fstr.append(f"{tbl} ")
        
        fstr = "".join(fstr)
        print(f"WRITING OUTPUT {output}")
        os.system(f"{fstr} > {output}")

    else:
        for tbl in tbls:
            df = pd.read_csv(tbl, sep='\t')
            df['batch'] = re.findall(r'batch[0-9]+', tbl)[0].split('batch')[1]
            dfList.append(df)

        combo = pd.concat(dfList)
        print(f"WRITING OUTPUT {output}")
        combo.sort_values(by = ['ParamGroup']).to_csv(output, sep = '\t', index=False)
        for tbl in tbls:
            os.rename(tbl, tbl.replace("code/CuBids","code/CuBIDS/batch"))


def combine(baseDir):
    
    summaries = []
    files = []
    scripts = []
    cubidsDir = op.join(baseDir, 'BIDS', 'code', 'CuBIDS')
    
    for fn in glob.glob(f"{cubidsDir}/*"):
        
        if '.sh' in fn:
            scripts.append(fn)
        elif re.search(r'run02_batch[0-9]+_summary.tsv', fn):
            summaries.append(fn)
        elif re.search(r'run02_batch[0-9]+_files.tsv', fn):
            files.append(fn)

    print(f"""COMBINING {len(scripts)} CMD SCRIPTS""")
    util = "ls "+cubidsDir+"/run02_batch*sh | xargs -I{} sed 1d {} > "+ op.join(cubidsDir, 'run02_full_cmd.sh')
    os.system(util)
    
    print(f"""COMBINING {len(summaries)} SUMMARY TABLES""")
    concat_tables(summaries, op.join(cubidsDir, 'run02_summary.tsv'))
    
    print(f"""COMBINING {len(files)} FILES TABLES""")
    concat_tables(files, op.join(cubidsDir, 'run02_files.tsv'))


def relaunch(stage, wrkDir, config, bidsDir, partial):

    """
    Relaunch failed jobs
    """

    # Load base config
    configDf = pd.read_csv(config)

    # Load job summary
    if stage == 'bidsify':
        df = pd.read_csv(op.join(wrkDir, 'curation_job_summary.tsv'), sep = "\t")
    else:
        df = pd.read_csv(op.join(wrkDir, 'postCubids_job_summary.tsv'), sep = "\t")

    # Make a new config based on the failed jobs
    failed_subs = df[df['status']=='FAILED']['pat_id'].unique().tolist()
    
    filtered_configDf = configDf[configDf['pat_id'].isin(failed_subs)]
    filtered_config = op.join(wrkDir, 'relaunch_config.csv')
    filtered_configDf.to_csv(filtered_config, index=False)

    n_jobs = df[df['status']=='FAILED']['pat_id'].nunique()
    print(f"RELAUNCHING {n_jobs} FAILED JOBS")
    print()

    # move any previous output to another directory, won't delete automatically just in case
    safe = op.join(wrkDir, 'relaunched_jobs_previous_output')
    if not op.exists(safe):
        os.mkdir(safe)

    print(f"MOVING PREVIOUS OUTPUT TO {safe}")
    print()
    for sub in failed_subs:
        print(sub)
        bids_dir = glob.glob(f"{bidsDir}/sub-{sub}")
        print(glob.glob(f"{bidsDir}/sub-{sub}"))
        heudiconv_trackers = glob.glob(f"/mnt/isilon/bgdlab_processing/Data/tmp_heudiconv/{sub}*")
        path = bidsDir.replace("BIDS","problematic_data")
        problematic_data = glob.glob(f"{path}/sub-{sub}")
        path = bidsDir.replace("BIDS","sourcedata")
        src_data = glob.glob(f"{path}/{sub}")

        lst = [item for sublist in [bids_dir, heudiconv_trackers,problematic_data, src_data] for item in sublist]
        for i in lst:
            new_i = i.replace(op.dirname(i), safe)
            print(f"mv {i} {new_i}")
            os.rename(i, new_i)
            
        new_src = op.join(bidsDir.replace("BIDS","sourcedata"), sub)
        os.mkdir(new_src)
        os.system(f"chmod 777 -R {new_src}")
            
    
    result, n_jobs, jobid = data_curate(wrkDir, filtered_config, bidsDir, partial)

    

def main():

    '''
    The main function of the script
    '''

    args = cli().parse_args()
    wrkDir = args.work_dir
    cmd = args.cmd
    baseDir = args.base
    delivery = args.dataset_name
    batch = args.dicoms
    pattern = args.pattern
    fnMap = args.cubids_full_cmd

    deliveryPath = op.join(baseDir, delivery)
    srcDir = op.join(deliveryPath, "sourcedata")
    bidsDir = op.join(deliveryPath, 'BIDS')


    # Check job status: --------------------------------------------
    if args.job_summary is True:

        if cmd == 'bidsify':
            summary = op.join(wrkDir, 'curation_job_summary.tsv')
            get_job_summary(summary)
        elif cmd == 'post-cubids':
            summary = op.join(wrkDir, 'postCubids_job_summary.tsv')
            get_job_summary(summary)
        sys.exit()


    # Prep and subject level organization: --------------------------------------------
    elif cmd == 'setup' or cmd == 'bidsify':

        # Sanity checks and preparations 
        for b in batch:
            
            if op.exists(b):
                pass
            else:
                print(f"Input directory not found {b}")
                sys.exit(1)

        deliveryExists = False if op.exists(deliveryPath) is False else True
        config = op.join(wrkDir, 'base_config.csv') if args.test is False else op.join(deliveryPath, 'base_config.csv')
        


        # Prep data
        if not op.exists(config) or cmd == 'setup':

            # Make the delivery directory
            if deliveryExists is False: 
                os.mkdir(deliveryPath)
            else:
                input(f"{deliveryPath} exists. Before executing commands to bidsify data, make sure that there is nothing in this path that could be potentially overwrriten.Press enter to continue ")
    

            result_prep = data_setup(wrkDir, args.manifest, " ".join(batch), pattern, config)
            print_cmd(result_prep)

            if result_prep.returncode == 0 or args.test is True:
                print()
                if result_prep.returncode == 0:
                    print(result_prep.stdout)
                
                # Create sourcedata directory 
                
                src = dl.Dataset(srcDir)
                src.create(cfg_proc='text2git', description = "Edited dicoms")
                bids = dl.Dataset(bidsDir)
                bids.create(cfg_proc='text2git', description = "BIDS compliant clinical imaging data")
                    
                
                # Create subject directories
                for subj in pd.read_csv(config).pat_id.unique():

                    os.mkdir(op.join(srcDir, subj))

                
                if src.is_installed():
                    src.save(path=srcDir, message = "Added subject level directories")

            else:
                print('FAILED TO SETUP DATA FOR CURATION')
                print(result_prep.stderr)
                sys.exit(1)

        # Subject level organization: --------------------------------------------
        if cmd == 'bidsify':

            if args.relaunch is False:
                # Sort dataframe so the subjects with the most sessions get organized first
                configDf = pd.read_csv(config)
                configDf = configDf.sort_values(by='original_subject_id', key=lambda x: x.str.len(), ascending=False)
                configDf.to_csv(config, index=False)

                result, n_jobs, jobid = data_curate(wrkDir, config, bidsDir, args.partial)
            else:
                relaunch('bidsify', wrkDir, config, bidsDir, args.partial)

    # CuBIDS: Top level organization: --------------------------------------------
    elif cmd == 'cubids-prep': 
        
        # Clean up processing files
        aggDfFn = op.join(wrkDir, "proc_summary.tsv")
        
        if not op.exists(aggDfFn):
            print("CLEANING UP PROCESSING FILES FROM STEP 1")
            # Remove txt files created by heuristic file used by heudiconv 
            configDf = pd.DataFrame(pd.read_csv(op.join(wrkDir, 'base_config.csv')))
            subs = configDf.original_subject_id.unique()
            filesToDelete = [glob.glob('/mnt/isilon/bgdlab_processing/Data/tmp_heudiconv/'+i+"*txt")[0]
                        for i in subs if (glob.glob('/mnt/isilon/bgdlab_processing/Data/tmp_heudiconv/'+i+"*txt"))]
            
            for f in filesToDelete: 
                os.remove(f)
            
            # Aggregate output tsv files
            outputs = [op.join(wrkDir,"logs","heudiconv"), op.join(deliveryPath, 'problematic_scans')]
            dfList = []
            for output in outputs:
                dfList = []
                if op.exists(output):
                    tables = glob.glob(output+f"/*.tsv")
                    for table in tables:
                        df = pd.DataFrame(pd.read_csv(table, sep = '\t'))
                        if 'heudiconv' in output and 'final_run' not in df.columns:
                            df['final_run'] = df['first_run']
                        dfList.append(df)
                    aggDf = pd.concat(dfList) if len(dfList) > 1 else df
                    aggDfFn = aggDfFn if 'heudiconv' in output else op.join(wrkDir, "problematic_scans.tsv")
                    aggDf.to_csv(aggDfFn, sep = "\t", index = False)

                    if not op.exists(aggDfFn):
                        sys.exit(1)
                
                    ## Remove individual tables
                    os.system(f"rm {output}/*''.tsv")
        
        # CuBIDS prep
        cubids_prep(bidsDir)
        


    # CuBIDS apply: --------------------------------------------
    elif cmd == 'cubids-apply':
        cubids_apply(wrkDir, deliveryPath, bidsDir)


    # Post CuBIDS sorting: --------------------------------------------
    elif cmd == 'post-cubids': 
        config = op.join(wrkDir, 'base_config.csv')
        post_cubids(wrkDir, config, fnMap, bidsDir)
    

    # Finish curation: --------------------------------------------
    elif cmd == 'annotate':
        annotate(wrkDir, bidsDir, baseDir, delivery, batch, args.manifest)


    # Combine batched CuBIDS outputs: --------------------------------------------
    elif cmd == 'combine':
        combine(baseDir)


if __name__ == "__main__":
    main()