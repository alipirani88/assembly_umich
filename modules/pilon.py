__author__ = 'alipirani'

import os
import subprocess
import re
import os
import errno
import glob
from config_settings import ConfigSectionMap
from modules.log_modules import keep_logging
from modules.logging_subprocess import *


def pilon(out_sorted_bam, reference, out_path, analysis, files_to_delete, logger, Config):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("pilon", Config)[
        'pilon_bin'] + "/" + ConfigSectionMap("pilon", Config)['base_cmd']
    keep_logging('Running Pilon for post-assembly improvement', 'Running Pilon for post-assembly improvement', logger, 'info')
    cmd = "java -jar %s --genome %s --bam %s --output %s --outdir %s --changes --verbose" % (
    base_cmd, reference, out_sorted_bam, analysis, out_path)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error at Pilon step. Exiting.',
                     'Error at Pilon step. Exiting.', logger, 'exception')
        sys.exit(1)

    return "%s/%s.fasta" % (out_path, analysis)