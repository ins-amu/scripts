import os
import shutil
import subprocess
import errno
import random
import string

def randomword(length=10):
   letters = string.ascii_lowercase
   return ''.join(random.choice(letters) for i in range(length))


def source_config_file(path_test_subj, config_file_name):
    # TODO: use subprocess.run from python 3.5
    command = ['bash', '-c', 'source '+ path_test_subj + '/' + config_file_name +' && env']
    print(command)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    # proc = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
    for line in proc.stdout:
      (key, _, value) = line.decode().partition("=")
      os.environ[key] = value
    proc.communicate()
    return dict(os.environ)


def get_path_parameters(path_test_subj, config_file_name):
# get path parameters from config file
    os.environ = source_config_file(path_test_subj, config_file_name)
    PRD = os.environ['PRD']
    SUBJ_ID = os.environ['SUBJ_ID']
    MATLAB = os.environ['MATLAB']
    return (PRD, SUBJ_ID, MATLAB)


def copy_dir(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

def setup_function(function):
    pass


def test_1():
    path_test_subj = "/disk2/Work/Processed_data/brown/scripts/cf-p08"
    test_dir_name = randomword()
    config_file_name = 'config_cf-p08.sh'
    PRD, SUBJ_ID, MATLAB = get_path_parameters(path_test_subj, config_file_name)
    env_test = source_config_file('tests/test_config_files/', 'default_config_file.sh') 
    env_test['PRD'] = PRD[:-1] + 'test_' +test_dir_name
    print(env_test['PRD'])
    env_test['SUBJ_ID'] = SUBJ_ID
    env_test['MATLAB'] = MATLAB
    env_test['SUBJECTS_DIR'] = "/disk2/Work/Processed_data/these/freesurfer"
    env_test['FREESURFER_HOME'] = "/usr/local/freesurfer"
    os.makedirs(os.path.join(path_test_subj, 'test_' + test_dir_name))
    copy_dir(os.path.join(path_test_subj, 'data'), os.path.join(path_test_subj, 'test_' + test_dir_name, 'data'))
    command = ['bash', './main_surface.sh', '-c', 'test']
    #command = ['bash', 'main_surface.sh', '-c', '/disk2/Work/Processed_data/brown/scripts/cf-p08/config_cf-p08.sh']
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, env=env_test)
    #proc = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
    #for stdout_line in iter(proc.stdout.readline, ""):
    #    yield stdout_line 
    #proc.stdout.close()
    #return_code = proc.wait()
    #if return_code:
    #    raise subprocess.CalledProcessError(return_code, cmd)
    for line in proc.stdout:
        print(line)
    proc.communicate()
    #assert proc.returncode==0

test_1()