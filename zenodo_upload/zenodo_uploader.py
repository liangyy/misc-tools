import yaml
import requests

def load_token(filename):
    with open(filename, 'r') as f:
        for i in f:
            i = i.strip()
            break
    return i
def load_yaml(ff):
    with open(ff, 'r') as f:
        e = yaml.safe_load(f)
    return e

def check_nfile(depository, token):
    token = load_token(token)
    params = {'access_token': token}
    headers = {"Content-Type": "application/json"}
    r = requests.get(
        'https://zenodo.org/api/deposit/depositions/{}'.format(depository), 
        params=params, json={}, headers=headers
    )
    tmp = r.json()
    if 'files' not in tmp:
        print('No file')
    else:
        print('{} file(s)'.format(len(tmp['files'])))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='zenodo_uploader.py', description='''
        Upload a list of files to zenodo.
        It will create one zenodo deposite and upload all files there.
    ''')
    parser.add_argument('--depository', default=None, help='''
        If want to append to existing depository, set the ID here.
        Otherwise, new depository will be generated.
    ''')
    parser.add_argument('--meta-yaml', help='''
        A yaml file containing some meta data to fill along with the files.
    ''')
    parser.add_argument('--file-list', help='''
        The list of files to upload.
    ''')
    parser.add_argument('--token', help='''
        Token saved in a file. 
        Note: make sure to delete the token after the uploading is done.
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    import ntpath
    import json
    
    # steps:
    # 1. create a depository
    # 2. update the meta data
    # 3. extract the bucket link of the depository
    # 4. upload files to the bucket
    # skip step 1 and 2 if args.depository is not None
   
     
    token = load_token(args.token)
    params = {'access_token': token}
    headers = {"Content-Type": "application/json"}
    if args.depository is None:
        # step 1
        r = requests.post("https://zenodo.org/api/deposit/depositions", params=params, headers=headers, json={})
        latest_draft = r.json()['links']['latest_draft']
        bucket_url = r.json()['links']['bucket']
        
        # step 2
        meta = load_yaml(args.meta_yaml)
        r = requests.put(latest_draft, params=params, data=json.dumps(meta), headers=headers)
    else:
        latest_draft = "https://zenodo.org/api/deposit/depositions/{}".format(args.depository)
        r = requests.get(latest_draft, params=params, headers=headers, json={})
        bucket_url = r.json()['links']['bucket']
        
    logging.info('Use depository = {}'.format(latest_draft))
    
    # step 3
    with open(args.file_list, 'r') as f:
        for i in f:
            i = i.strip()
            filename = ntpath.basename(i)
            logging.info(f'Processing {filename}')
            with open(i, "rb") as fp:
                r = requests.put(
                    "%s/%s" % (bucket_url, filename),
                    data=fp,
                    params=params,
                )
            
