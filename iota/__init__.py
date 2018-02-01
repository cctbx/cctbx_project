from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 01/29/2018
Description : IOTA initialization module (also contains app info)
'''

from datetime import datetime


iota_version = '1.1.026'
now = "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())

# For GUI
gui_description = '''The integration optimization, triage and
analysis (IOTA) toolkit for the processing of serial diffraction data.

Reference: Lyubimov, et al., J Appl Cryst, 2016
'''

prime_description = ''' The Post-RefInement and MErging (PRIME) program for
the scaling, merging and post-refinement of integrated diffraction images.

Reference: Uervirojnangkoorn, et al., eLife, 2015'''

gui_license = ''' IOTA is distributed under open source license '''
prime_license = ''' PRIME is distributed under open source license'''
cxi_merge_license = ''' cxi.merge is distributed under open source license '''
