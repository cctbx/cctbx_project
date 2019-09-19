from __future__ import absolute_import, division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 09/17/2019
Description : IOTA initialization module (also contains app info)
'''

from datetime import datetime

iota_version = '1.4.009'
intx_version = '1.0.001'
now = "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())

# For GUI
gui_description = '''The integration optimization, triage and
analysis (IOTA) toolkit for the processing of serial diffraction data.

Reference: Lyubimov, et al., J Appl Cryst, 2016
'''

gui_license = ''' IOTA is distributed under open source license '''
cxi_merge_license = ''' cxi.merge is distributed under open source license '''

help_message = '\n{:-^70}'\
               ''.format('Integration Optimization, Triage and Analysis') + """

DRY RUN - command line
Usage: iota.run -d
Will output an IOTA script file (iota.param) and a cctbx.xfel
target file (cctbx_xfel.phil), and also show the IOTA script
with help strings in terminal window.

COMMAND LINE MODE
Usage: iota.run [FILES] [OPTIONS] [PHIL ARGS]
IOTA accepts a variety of input types, including paths to raw
images, individual imagefiles, IOTA script file, text files listing
imagefile paths, or any combination thereof. If both an IOTA
script file and paths to images are specified, the image lists
will be combined and purged of duplicates. IOTA settings can be
altered using command line arguments (see examples below). If a
script file is not provided, defaults will be generated automatically.

Examples:

iota.run /Users/doug/cool_hewl_data/ -n 10 image_import.mask=mask.pickle

iota.run iota.param ../cool_hewl_data/*.cbf -n 10 image_import.mask=None

iota.run iota.param -n 10 image_import.mask=mask.pickle


GUI MODE
Usage: iota [FILES] [OPTIONS] [PHIL ARGS]
This will launch the IOTA GUI; the GUI will be automatically
populated with file paths and settings if they are provided.
Otherwise, a default configuration with no data will be used.

Examples:

iota                             # will launch the GUI with defaults

iota /Users/doug/cool_hewl_data/ -n 10 image_import.mask=mask.pickle

iota iota.param ../cool_hewl_data/*.cbf -n 10 image_import.mask=None

"""

logo = """
     IIIIII            OOOOOOO        TTTTTTTTTT          A
       II             O       O           TT             A A
       II             O       O           TT            A   A
>------INTEGRATION----OPTIMIZATION--------TRIAGE-------ANALYSIS------------>
       II             O       O           TT          A       A
       II             O       O           TT         A         A
     IIIIII            OOOOOOO            TT        A           A   v{}
"""
logo = logo.format(iota_version)

# -- end
