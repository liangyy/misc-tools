# upload file to Box

import argparse
parser = argparse.ArgumentParser(prog='upload_to_box.py', description='''
    Upload file to Box
''')
parser.add_argument('--client_id', required=True, help='''
    Client ID for boxsdk.
''')
parser.add_argument('--client_secret', required=True, help='''
    Client secret for boxsdk.
''')
parser.add_argument('--developer_token', required=True, help='''
    Developer token for boxsdk.
''')
parser.add_argument('--folder_id', required=True, help='''
    Folder ID of the target Box folder.
''')
parser.add_argument('--file_path', required=True, help='''
    Path to file to be uploaded.
''')
args = parser.parse_args()

from boxsdk import Client, OAuth2

oauth = OAuth2(
    client_id=args.client_id,
    client_secret=args.client_secret,
    access_token=args.developer_token,
)

client = Client(oauth)

folder_id = args.folder_id

print('Start uploading {}'.format(args.file_path))
new_file = client.folder(folder_id).upload(args.file_path)
print('Finish uploading {}'.format(args.file_path))

