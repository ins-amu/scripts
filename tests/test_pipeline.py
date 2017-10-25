import os
import shutil
import subprocess


def source_config_file(path, config_file_name):
# TODO: use subprocess.run from python 3.5
    command = ['bash', '-c', 'source '+ path_test_subj + ' ' + config_file_name +' && env']
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    # proc = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
    for line in proc.stdout:
      (key, _, value) = line.decode().partition("=")
      print(value)
      print(key)
      os.environ[key] = value
    proc.communicate()
    return dict(os.environ)
 

def get_path_parameters(path_test_subj, config_file_name):
# get path parameters from config file
    os.environ = source_config_file(path_test_subj, config_file_name)
    PRD = os.environ['PRD']
    SUBJ_ID = os.environ['SUBJ_ID']
    MCR = os.environ['MCR']
    return (PRD, SUBJ_ID, MCR)


def copy_dir(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

def setup_function(function):
    path_test_subj = ""
    config_file_name = 'config_ac-p15.sh'
    PRD, SUBJ_ID, MCR = get_path_parameters(path_test_subj, config_file_name)


def test_1(self):
    os.environ = source_config_file('./test_config_files/', 'default_config_file.sh')  
    os.environ['PRD'] = PRD
    os.environ['SUBJ_ID'] = SUBJ_ID
    os.environ['MCR'] = MCR
    os.makedirs(os.path.join(path_test_subj, 'test_1'))
    copy_dir(os.path.join(path_test_subj, 'data'), os.path.join(path_test_subj, 'test_1')
    command = ['bash', 'main_surface.sh', '-c', 'test']
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    for line in proc.stdout:
        print(line)
    assert proc.returncode==0
