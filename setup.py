import sys
import logging
import subprocess
import os
import tempfile
import shutil
import logging
try:
    import numpy
except ModuleNotFoundError:
    m = "ERROR. Please install numpy in your Python environment.\n"
    m += "ERROR. pip install numpy"
    print(m)
    exit()
from datetime import datetime
from setuptools import setup, Extension
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler


# Formatter for the logger
class CustomFormatter(logging.Formatter):

    """Logging Formatter to add colors and count warning / errors"""
    FORMATS = {
        logging.ERROR: "\n\tERROR: %(asctime)s: %(msg)s",
        logging.WARNING: "\n\tWARNING: %(msg)s",
        logging.DEBUG: "%(asctime)s: %(msg)s",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        date_fmt = '%d-%m-%Y %d %H:%M:%S'
        formatter = logging.Formatter(log_fmt, date_fmt)
        return formatter.format(record)


# Install packages from pip ==============================================================
def install_with_pip(pack, vers=None, log=None, namepkg=None):

    # sys.executable gives the path of the python interpreter
    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if vers is None:
        m = "{}: ** {}: Installing {}".format(namepkg, now, pack)
        print(m) if log is None else log.info(m)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}".format(pack)])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}".format(pack)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()
    else:
        m = "{}: ** {}: Installing {}=={}".format(namepkg, now, pack, vers)
        print(m) if log is None else log.info(m)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers), " &>install.log"])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()


# Main setup
if __name__ == '__main__':

    namepackage = "pyCGMapPol"

    # Creating the logger to install.log file ===================================
    logger = logging.getLogger(name="INSTALL_LOG")
    logger.setLevel(logging.DEBUG)
    h1 = logging.FileHandler("install.log", 'w')
    h1.setFormatter(CustomFormatter())
    # Output also in the screen
    logger.addHandler(h1)
    f1 = logging.StreamHandler()
    f1.setFormatter(CustomFormatter())
    logger.addHandler(f1)

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Starting installation!!!! at {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)

    # Print sys path ===================================
    m1 = "\t\t SYS PATH\n"
    for item in sys.path:
        m1 += item + "\n"
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 += "\n\t\t INSTALLING PIP PACKAGES ({})\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)

    # Install requirements ===================================
    with open('requirements.txt') as f:
        required = f.read().splitlines()
    for ipack in required:
        try:
            pkg, version = ipack.split(">=")[0:2]
            if pkg[0] == "#":
                continue
            install_with_pip(pkg, vers=version, log=logger, namepkg=namepackage)
        except ValueError:
            pkg = ipack
            try:
                if pkg[0] == "#":
                    continue
                install_with_pip(pkg, log=logger, namepkg=namepackage)
            except IndexError:
                pass

    m1 = "\n\t\t RUNNING SETUP FROM SETUPTOOLS {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    setup()

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Installation Done!!!! at {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)

