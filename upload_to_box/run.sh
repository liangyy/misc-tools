# upload models to Box

clientId=''
clientSecret=''
accessToken=''

uploadScript=upload_to_box.py

files='file1 file2'
fid=''  # folder id in Box

for f in $files
do
  python $uploadScript \
  --client_id $clientId \
  --client_secret $clientSecret \
  --developer_token $accessToken \
  --folder_id $fid \
  --file_path $f
done
