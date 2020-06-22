#!/bin/bash
#
# Apply transformations from atlas image to the anatomical image
# Fabien Boux 04/2018

# ----- Arguments -----
DIM=3
VERBOSE=0

while getopts :s:l:d:t:i:o:v option
do
  case ${option}
    in
    s) SCRIPT_LOC="$OPTARG"
        ANTSPATH="$OPTARG"
        ;;
    l) LABEL_IMG="$OPTARG"
        ;;
    d) DIM="$OPTARG"
        ;;
    t) TRANSFO="$OPTARG"
        ;;
    i) IMAGE="$OPTARG"
        ;;
    o) OUTPUT="$OPTARG"
        ;;
    v) VERBOSE="$OPTARG"
        ;;
  esac
done

if [ "$SCRIPT_LOC" == "" ] || [ "$LABEL_IMG" == "" ] || [ "$IMAGE" == "" ];
  then
    echo "Error : Something wrong with arguments"
    exit 2
fi


export ANTSPATH=$ANTSPATH


# ----- Script -----

# Apply transformations
echo -e "------------------------------------------------------------------------\n\
TODO \n\
 Image : $IMAGE \n\
 Dimension: $DIM \n\
------------------------------------------------------------------------"

$SCRIPT_LOC/antsApplyTransforms \
                                 -d $DIM \
                                 -t  $TRANSFO \
                                 -r  $IMAGE \
                                 -i  $LABEL_IMG \
                                 -o $OUTPUT \
                                 -n NearestNeighbor \
                                 -v $VERBOSE
