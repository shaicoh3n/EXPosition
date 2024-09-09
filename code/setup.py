import sys
from cx_Freeze import setup, Executable
import os

include_files = [os.path.join('.','EXPosition_utils.py'),
                 os.path.join('.','Data'),
                 os.path.join('.','Models'),
                 os.path.join('.','transcription'),
                 os.path.join('.','annotations'),
                 os.path.join('.','Results'),
                 os.path.join('.','man_muts_example.csv'),
                 os.path.join('.','pred_muts_example.csv')]

packages = ['numpy', 'pandas','sklearn','scipy','tensorflow','keras','statsmodels','tqdm','termplotlib']
setup(
    name='EXPostion',
    version='1.0',
    options={
        'build_exe': {
            'include_files': include_files,
            'packages': packages,
        }
    },
    executables=[Executable('Exposition_GUI.py')]
)