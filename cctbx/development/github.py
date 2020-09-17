from __future__ import absolute_import, division, print_function

import json
import os
import requests
import time

from libtbx.program_template import ProgramTemplate

class Program(ProgramTemplate):

  description = '''
Create an issue on GitHub via their REST API.
'''

  master_phil_str = '''
  org = None
    .type = str
    .help = GitHub organization
    .short_caption = GitHub organization
  repo = None
    .type = str
    .help = GitHub repository
    .short_caption = GitHub repository
'''

def create_github_issue(org=None, repo=None, title='title', body='body', token=None):
  assert token is not None
  url = 'https://api.github.com/repos/{org}/{repo}/issues'.format(org=org, repo=repo)

  # headers with authenticating token
  headers = {
    'Accept': 'application/vnd.github.v3+json',
    'Authorization': 'token {token}'.format(token=token),
  }

  # issue
  issue = {
            'title': title,
            'body': body,
          }
  data = json.dumps(issue)

  # create issue
  response = requests.post(url, data=data, headers=headers)
  if response.status_code == 201:
    print('Successfully created issue')
  else:
    print('Issue was not created')

  return response

if __name__ == '__main__':
  org = 'phenix-lbl'
  repo = 'testing'
  title = 'test {}'.format(time.asctime())
  body = 'This is a test {}'.format(time.asctime())
  home = os.path.expanduser('~')
  with open('{home}/github.token'.format(home=home), 'r') as f:
    token = f.read()
  response = create_github_issue(org=org, repo=repo, title=title, body=body, token=token)
