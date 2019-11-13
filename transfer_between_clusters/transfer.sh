# ARGV1: local path to folder that you would like to transfer (folder with / at the end)
# ARGV2: the target server user@server:path
# ARGV3: log file name
# What it does: copy folder to server path (it should be a NEW folder otherwise it will put the content in existing folder)
# Before run: you should make sure that you have access to server (e.g. via ssh key, etc)

infolder=$1
case $infolder in      */) echo The input folder ends with / ;;      *) infolder=$infolder/ ;; esac

rsync -avc --log-file=$3 $infolder $2