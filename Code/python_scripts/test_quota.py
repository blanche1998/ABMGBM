import subprocess
import numpy as np
import os
import sys
import shutil #this will be used to copy files by using the command shutil.copyfile(src, dst)
              #in that case we are copying the content of the file src into the file dst
import xml.etree.ElementTree as ET #to parse through and modify a xml file
from datetime import datetime

#NOBACKUP_path = "/NOBACKUP/mongeonb/"

#total, used, free = shutil.disk_usage(NOBACKUP_path)

#print("Total: %d GiB" % (total // (2**30)))
#print("Used: %d GiB" % (used // (2**30)))
#print("Free: %d GiB" % (free // (2**30)))

cmd = 'quota'

os.system(cmd)