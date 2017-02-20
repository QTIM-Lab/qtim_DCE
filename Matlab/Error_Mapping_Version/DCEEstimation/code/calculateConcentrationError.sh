#!/bin/bash


#-----------------------------
help()
{
cat << HELP

This script calculate the concentration estimation errors

DEPENDENCIES:

	fslsplit	: `which fslsplit`
	fslmerge	: `which fslmerge`
	fslmaths	: `which fslmaths`


########################################
USAGE: $0 [OPTIONS]
OPTIONS:

Reqd:   -in1   <file>  : observed  concentration (4D file)
	-in2   <file>  : estimated concentration (4D file)

Optional
	-roi   <file>   roi file (to calculate difference only within roi)
	

HELP
exit 1
}


#-----------------------------
parse()
{
	while [ -n "$1" ]; do
		case $1 in
			-h) 
				help;
				shift 1;;			# help is called
			
		     	-in1) 
				input1=$2;
				shift 2;;			# input image 1

			-in2) 
				input2=$2;
				shift 2;;			# input image 2

			-roi)
				roi=$2;
				shift 2;;

			-*) 
				echo "ERROR: no such option $1";
				helpshort;;
			 *) 
				break;;
		esac
	done
}



#---------------------------------
# load library or executables, they are dependencies of this pipeline
loaddependency()
{
	AFNI_DIR=/usr/pubsw/packages/AFNI/current
	FSL_DIR=/usr/pubsw/packages/fsl/current/bin
	for dir in ${AFNI_DIR} ${FSL_DIR}
	do
	    if [ -d ${dir} ]; then
	      export PATH=$PATH:${dir} 
	    fi
	done
}


if [ $# -lt 2 ]
then
	help
fi


### Reading the arguments
parse $*


##  main 
echo ""
echo "calculate the DCE estimation differences between"
echo "          the observed  concentration ${input1}, and"
echo "          the estimated concentration ${input2}"
echo ""

diffsqrdimage4D=diffsqrt4D.nii.gz
meandiffimage3D=meandiff.nii.gz

# calculate the difference
fslmaths ${input1} -sub ${input2} -abs -sqrt ${diffsqrdimage4D}
fslmaths ${diffsqrdimage4D} -Tmean ${meandiffimage3D}

if [ ! -z ${roi} ]; then
    dmean=`3dmaskave -mask ${roi}       ${meandiffimage3D} |awk '{ print $1} '`  # mean squared differences
    dmin=` 3dmaskave -mask ${roi} -min  ${meandiffimage3D} |awk '{ print $1} '`  # min  squared differences
    dmax=` 3dmaskave -mask ${roi} -max  ${meandiffimage3D} |awk '{ print $1} '`  # max  squared differences
else
    dmean=`3dmaskave       ${meandiffimage3D} |awk '{ print $1} '`  # mean squared differences
    dmin=` 3dmaskave -min  ${meandiffimage3D} |awk '{ print $1} '`  # min  squared differences
    dmax=` 3dmaskave -max  ${meandiffimage3D} |awk '{ print $1} '`  # max  squared differences
fi

\rm ${diffsqrdimage4D} ${meandiffimage3D}


# output
echo ""
echo "mean_squared_difference = ${dmean}"
echo "max__squared_difference = ${dmax}"
echo "min__squared_difference = ${dmin}"
echo ""

