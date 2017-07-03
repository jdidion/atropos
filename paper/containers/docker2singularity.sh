#!/bin/bash
# Usage: docker2singularity.sh <docker2singularity container name> <remote host> <remote host dir>
D2S=$1
REMOTE_HOST=$2
REMOTE_DIR=$3
while read -r container
do
   echo "Transferring $container to ${REMOTE_HOST}:${REMOTE_DIR}" \
&& name=`basename $container` \
&& file="${name}.img" \
&& docker run \
   -v /var/run/docker.sock:/var/run/docker.sock \
   -v $(pwd):/output --privileged -t --rm $D2S $container \
&& mv *${name}*.img $file \
&& echo "scp $file ${REMOTE_HOST}:${REMOTE_DIR}" \
&& scp $file ${REMOTE_HOST}:${REMOTE_DIR} 
done < <(grep -v -P "^#" paper_containers.txt)
