#!/bin/bash
# Starts a NIAK instance inside a singularity container 
# and a process to execute code on the host, typically qsub


CONFIG_FILE=psom.conf
#Options for singularity execution
PSOM_SINGULARITY_OPTIONS=''

# Were singularity images are keept 
PSOM_SINGULARITY_IMAGES_PATH=${HOME}/singularity/:.

# Directory contaning .m files to be loaded at startup 
#    typically contains "psom_gb_vars_local.m" 
PSOM_LOCAL_CONF_DIR='' 

# load config file
# Gobal
source /etc/${CONFIG_FILE} > /dev/null 2>&1
# Global when no su acces in granted
CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${CURRENT_DIR}/${CONFIG_FILE}
#User specific
source ${HOME}/.config/psom/${CONFIG_FILE} > /dev/null 2>&1


list_all_image () {

    while IFS= read -r -d $'\0' line; do
        if head -n1 $line | grep -q run-singularity ; then
          ALL_IMAGES=${line}:${ALL_IMAGES%:}
        fi
    done  < <(find ${PSOM_SINGULARITY_IMAGES_PATH//:/ } -maxdepth 1 -type f -print0)

    # Available images
    echo Installed images:
    while read -r line ; do
      ONE_IMAGE=
      bidon=${line##*/}
      echo "     ${bidon%.img}"
    done < <(echo ${ALL_IMAGES}| tr ':' '\n')
    if [[ ! ${ONE_IMAGE+x} ]]; then
       echo No Image installed
       echo "try running with the -p <path_to_image> option"
    fi
}

usage (){

    echo "Starts octave in psom/singulrity mode"
    echo
    echo "Usage: $(basename $0)  <installed_image>"
    echo "   or: $(basename $0)  -p <path_to_singularity_image>"
    echo "   or: $(basename $0)  -l"
    echo
    echo "   -l                  List locally installed images"
    echo
    echo "   -p                  Run from an arbitrary singularity image path"
    echo
}

finish () {
  # Your cleanup code here
  kill ${LOOP_ID} > /dev/null 2>&1
  rm -r ${PSOM_FIFO_DIR} > /dev/null 2>&1
  rm -r /tmp/psom_worker_qsub.* > /dev/null 2>&1
  ## Could add a qdel cmd here to clean process when console is not running
}
trap finish EXIT

host_exec_loop () {
# creates a tmp file to execute code on singularity via qsub
  while true
  do
    if read line; then
       QSUB_OPT=(${line%SPLIT_LINE*})
       if [[ ${QSUB_OPT[0]} == 'qsub_options' ]]; then
         SINGULARITY=(${line#*SPLIT_LINE})
         tmpfile=$(mktemp /tmp/psom_worker_qsub.XXXXXX)
         echo '#!/bin/bash' > ${tmpfile}
         echo singularity  exec ${PSOM_SINGULARITY_OPTIONS} ${SINGULARITY[@]:1}  >> ${tmpfile}
         qsub ${QSUB_OPT[@]:1} ${tmpfile}
       else
           echo INVALID qsub CMD "${QSUB_OPT[@]}"
     fi
   fi
  done <"$PSOM_FIFO"
}


# Check opt args
while getopts ":lp:" opt; do
  case $opt in
    l)
      list_all_image
      exit 0
      ;;
    p)
      IMAGE_PATH=$OPTARG
      if [ ! -f "${IMAGE_PATH}" ]; then
        echo ${IMAGE_PATH} not found
        usage
      fi
      ;;
    \?)
      usage
      list_all_image
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      list_all_image
      exit 1
      ;;
  esac
done

# Check pos args
if [ -z ${IMAGE_PATH} ] ; then
    shift $(expr $OPTIND - 1 )

    while test $# -gt 0; do
      pos_arg=$1
          IMAGE_PATH=$(find ${PSOM_SINGULARITY_IMAGES_PATH//:/ } \
                        -maxdepth 1 -type f -name "${pos_arg%.img}.img" -print -quit)
          if [[ -z "${IMAGE_PATH// }" ]]; then
            echo image ${OPTARG} not found
            echo
            echo Please select
            list_all_image
            echo
            usage
            exit 1
          else
              break
          fi
    done
fi

if [ -z ${IMAGE_PATH} ] ; then
  usage
  list_all_image
  exit
fi

IMAGE_PATH=$(readlink -f ${IMAGE_PATH})

CONSOLE_ID=$$

export PSOM_FIFO_DIR=$(mktemp -d /tmp/psom-${CONSOLE_ID}-fifo.XXXXXX)
export PSOM_FIFO=${PSOM_FIFO_DIR}/pipe

[ -p $PSOM_FIFO ] || mkfifo $PSOM_FIFO;

# Start the communication loop
host_exec_loop > /dev/null 2>&1 & 
LOOP_ID=$!

# Start singularity-psom

PSOM_START_OCTAVE="octave --persist --no-init-file --eval \"addpath(genpath(\\\"${PSOM_LOCAL_CONF_DIR}\\\"))\"" 
singularity shell ${PSOM_SINGULARITY_OPTIONS} ${IMAGE_PATH} -c "export PSOM_FIFO=${PSOM_FIFO};export PSOM_LOCAL_CONF_DIR=${PSOM_LOCAL_CONF_DIR}  ; ${PSOM_START_OCTAVE}"
