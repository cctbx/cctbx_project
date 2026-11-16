from __future__ import absolute_import, division, print_function
import wx

# Location of the GUI documentation on the cctbx.xfel MediaWiki. HELP_URL is the
# main page; each tab of the main window has its own subpage beneath it, named
# in HELP_TAB_PAGES, and dialogs name their subpages directly. Until HELP_URL is
# set to None, the Help buttons fall back to the published GUI description.
HELP_URL = 'https://cci.lbl.gov/xfel/index.php/Cctbx.xfel_GUI'
HELP_FALLBACK_URL = 'http://cci.lbl.gov/publications/download/CCN_2019_p22_Brewster.pdf'
HELP_TAB_PAGES = {
  'Runs':          'Runs_tab',
  'Energy':        'Energy_tab',
  'Trials':        'Trials_tab',
  'Jobs':          'Jobs_tab',
  'Run Stats':     'Run_Stats_tab',
  'Unit Cell':     'Unit_Cell_tab',
  'Datasets':      'Datasets_tab',
  'Merging stats': 'Merging_Stats_tab',
}

def help_url(page=None):
  ''' The documentation page named (a subpage of HELP_URL), or the main page
      when no page is named. '''
  if HELP_URL is None:
    return HELP_FALLBACK_URL
  if page is None:
    return HELP_URL
  return '%s/%s' % (HELP_URL.rstrip('/'), page)

def open_help(parent, page=None):
  ''' Open the named documentation page in the default web browser, showing
      the address instead when no browser could be opened. '''
  import webbrowser
  url = help_url(page)
  print('Opening', url)
  if not webbrowser.open(url):
    wx.MessageBox('Could not open a web browser. The documentation is at:\n%s' % url,
                  'Help', wx.OK | wx.ICON_INFORMATION, parent)
