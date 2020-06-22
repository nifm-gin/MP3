#!/bin/bash
#
# Script to match images with atlas
#
# Fabien Boux 03/2018


# ----- Arguments -----
DIM=3
TRANSFO=s

while getopts :s:a:i:d:t:m:o: option
do
  case ${option}
    in
    s) SCRIPT_LOC="$OPTARG"
        ANTSPATH="$OPTARG"
        ;;
    a) ATLAS="$OPTARG"
        ;;
    i) IMAGE="$OPTARG"
        ;;
    d) DIM="$OPTARG"
        ;;
    t) TRANSFO="$OPTARG"
        ;;
    m) MASK="$OPTARG"
        ;;
    o) OUTPUT="$OPTARG"
        ;;
  esac
done

VERBOSE=1

if [ "$SCRIPT_LOC" == "" ] || [ "$ATLAS" == "" ] || [ "$IMAGE" == "" ];
  then
    echo "Error : Something wrong with arguments"
    exit 2
fi


export ANTSPATH=$ANTSPATH


# Execute inhomogeneity field correction
#echo -e "------------------------------------------------------------------------\n\
#Apply inhomogeneity field correction \n\
# Image: $IMAGE \n\
# Mask: $MASK \n\
# Dimension: $DIM \n\
#------------------------------------------------------------------------"

#CORRECTED_IMG=$OUTPUT-corrected.nii.gz

#$SCRIPT_LOC/N4BiasFieldCorrection \
#                                 -d $DIM \
#                                 -i  $IMAGE \
#                                 -x $MASK \
#                                 -o $CORRECTED_IMG \
#                                 -v $VERBOSE


# Compute transformation
echo -e "------------------------------------------------------------------------\n\
Compute transformation \n\
 Image: $IMAGE \n\
 Atlas: $ATLAS \n\
 Dimension: $DIM \n\
 Output: $OUTPUT \n\
------------------------------------------------------------------------"

$SCRIPT_LOC/antsRegistrationSyNQuick.sh \
                                 -d $DIM \
                                 -f $ATLAS \
                                 -m $IMAGE  \
                                 -t $TRANSFO \
                                 -r  64 \
                                 -j 1 \
                                 -o $OUTPUT



