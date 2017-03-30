#!/bin/bash

upload_src_dir=$1
remote_dir=$2
admin_file=$3

source $admin_file
info=$upload_src_dir/Info

upload_sh=`pwd`/sendme.sh

#adding author's name
me=`whoami`
echo "Author:$me" >$info
echo "Date: `date`">>$info

cmd="`which upload.sh` `readlink -f $upload_src_dir` $remote_dir $remote_admin"

if [ `whoami` = $local_admin ]; then
    $cmd
else
    echo "!!!!!  Ask Admin to upload: Give him the path: `pwd` !!!!!"
    echo "#!/bin/bash" >$upload_sh
    echo $cmd >>$upload_sh
    chmod +x $upload_sh
fi


