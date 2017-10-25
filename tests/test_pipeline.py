import os
import shutil
import pprint
import subprocess


command = ['bash', '-c', 'source config_ac-p15.sh && env']

proc = subprocess.Popen(command, stdout = subprocess.PIPE)


for line in proc.stdout:
  (key, _, value) = line.decode().partition("=")
  print(value)
  print(key)
  os.environ[key] = value

proc.communicate()
pprint.pprint(dict(os.environ))

def copy_dir(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

def setup_function(function):
    path_test_subj = ""
    os.makedirs(os.path.join(path_test_subj, 'test_1'))
    copy_dir(os.path.join(path_test_subj, 'data'), os.path.join(path_test_subj, 'test_1')
    shutil.copyfile(exam


def test_1(self):





