#!/sw/bin/python2.6

"""
A sample program that runs whisker.py as an embedded function rather than 
a standalone executable.
"""

import logging
import sys

logging_level = logging.WARNING
#logging_level = logging.INFO
#logging_level = logging.DEBUG
logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
logger = logging.getLogger("whisker")

kwargs = {
    'root_dir' : '../../DECam_159069',
    'exp_num' : 159069,
    # The rest are optional:
    'out_file' : 'wl.dat',
    #'make_plots' : True,
    'logger' : logger
}
import whisker 
table = whisker.process_all(**kwargs)

print 'Returned table = ',table
