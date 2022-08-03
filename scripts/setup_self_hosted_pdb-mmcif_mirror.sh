#!/bin/bash
# edited by Yinying Yao
set -e
db_path=$(readlink -f $1)

echo "----------------------------------------------------------------------------"
echo "Checking for linux distribution ..."
if [[ "$(which yum)" != "" ]];then
  LINUX_DIST="REDHAT"
elif [[ "$(which apt)" != "" ]];then
  LINUX_DIST="DEBIAN"
else
  LINUX_DIST="DEBIAN"
fi
echo "----------------------------------------------------------------------------"
echo "Checking for rsync installation ..."
if [[ "$(which rsync)" == "" ]];then
  echo "Rsync not available."
  if [[ "${LINUX_DIST}" == "REDHAT" ]];then
    sudo yum install rsync -y
  elif [[ "${LINUX_DIST}" == "DEBIAN" ]];then
    sudo apt install rsync -y
  else
    echo "Error: Please install rsync by yourself."
    exit 1
  fi
fi

mkdir -p $db_path

rsync_configuration="# /etc/rsyncd: configuration file for rsync daemon mode \n

# See rsyncd.conf man page for more options.\n

uid = nobody\n
gid = nobody\n
use chroot = no\n
max connections = 0\n
lock file=/var/run/rsyncd.lock\n
log file = /var/log/rsyncd.log\n
exclude = lost+found/\n
transfer logging = yes\n
timeout = 900\n
ignore nonreadable = yes\n
dont compress   = *.gz *.tgz *.zip *.z *.Z *.rpm *.deb *.bz2\n
\n
[db]\n
path = ${db_path}\n
comment=dbs\n
ignore errors\n
read only=yes\n
write only=no\n
list=no\n
#auth users=user\n
#secrets file=/etc/rsyncd.passwd\n
\n
hosts allow=* \n"

echo "----------------------------------------------------------------------------"
echo "Now configure rsyncd ..."
echo "Your configurations are listed below:"
echo "----------------------------------------------------------------------------"
echo -e ${rsync_configuration}
echo "----------------------------------------------------------------------------"
sudo systemctl start rsyncd.service
sudo systemctl enable rsyncd.service
echo -e ${rsync_configuration} >/etc/rsyncd.conf
rsync --daemon --config=/etc/rsyncd.conf
sudo systemctl restart rsyncd.service

echo "----------------------------------------------------------------------------"
echo "Modifying firewall configuration in expected linux distribution ${LINUX_DIST} ... "
if [[ "$LINUX_DIST" == "CENTOS" ]]; then
  cmd="sudo firewall-cmd --zone=public --add-port=873/tcp --permanent"
  echo "$cmd" && eval "$cmd" || exit 1
elif [[ "$LINUX_DIST" == "DEBIAN" ]]; then
  cmd="sudo ufw allow 873 && sudo ufw reload"
  echo "$cmd" && eval "$cmd" || exit 1
fi

echo "----------------------------------------------------------------------------"
echo "Now we start to synchronize from original mirrors ..."
mkdir --parents "${db_path}/pdb_mmcif"

echo "Trying to sync with PDBj ..." && \
MIRROR="PDBj" && \
MIRROR_CMD="rsync --recursive --links --perms --times --compress --info=progress2 --delete data.pdbj.org::ftp_data/structures/divided/mmCIF/ ${db_path}/pdb_mmcif/raw/" && \
echo "${MIRROR_CMD}" && eval "${MIRROR_CMD}" || \
echo "Trying to sync with EBI ..." && MIRROR="EBI" && \
MIRROR_CMD="rsync --recursive --links --perms --times --compress --info=progress2 --delete rsync.ebi.ac.uk::ftp_data/structures/divided/mmCIF/ ${db_path}/pdb_mmcif/raw/" && \
echo "${MIRROR_CMD}" && eval "${MIRROR_CMD}" || \
echo "Trying to sync with RCSB PDB ..." && MIRROR="RCSB" && \
MIRROR_CMD="rsync --recursive --links --perms --times --compress --info=progress2 --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/mmCIF/  ${db_path}/pdb_mmcif/raw/" &&
echo "${MIRROR_CMD}" && eval "${MIRROR_CMD}" || echo "Failed to run all rsync tests" && exit 1

echo "----------------------------------------------------------------------------"
echo "Done! Mirror site: ${MIRROR}"
echo "The following command will be executed every Friday:"
echo "${MIRROR_CMD}"
echo "----------------------------------------------------------------------------"
echo "Now we create a crontab task."
CRONTAB_CMD="#!/bin/bash\n
PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/sbin:~/bin\n
export PATH\n
${MIRROR_CMD}\n
echo \"----------------------------------------------------------------------------\"\n
endDate=\`date +\"%Y-%m-%d %H:%M:%S\"\`\n
echo \"★[\${endDate}] Successful\"\n
echo \"----------------------------------------------------------------------------\"\n
"
echo "Your crontab task configuration is listed here:"
echo "----------------------------------------------------------------------------"
echo -e "${CRONTAB_CMD}"
echo "----------------------------------------------------------------------------"
echo "${CRONTAB_CMD}" > "${db_path}/pdb_mmcif_rsync.sh"
touch /var/spool/cron/$(whoami)
echo "30 16 * * 5 ${db_path}/pdb_mmcif_rsync.sh >> ${db_path}/pdb_mmcif_rsync.log 2>&1" >> /var/spool/cron/$(whoami)

echo "Congrats! Now you have a mirror site against the ${MIRROR}!"




