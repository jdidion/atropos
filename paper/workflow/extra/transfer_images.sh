REMOTE_HOST=$1
REMOTE_DIR="/scratch/jdidion/containers"
while read -r container
do
   echo "Transfering $container to ${REMOTE_HOST}:${REMOTE_DIR}" \
&& name=`basename $container` \
&& file="${name}.tar" \
&& docker save $container > $file \
&& scp $file ${REMOTE_HOST}:${REMOTE_DIR} \
&& rm -Rf $file
done < <(grep -v -P "^#" ../../containers/paper_containers.txt)
