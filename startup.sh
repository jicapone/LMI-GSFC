#
#
#!/bin/bash

##########################
# check for dependencies #
##########################
echo "Checking pipeline add-on dependencies:"

# check that pipeline has been sourced
[ -z "$PIPE_ROOT" ] && echo "Need to run the pipeline's startup.sh first.  Aborting." && exit 1
echo -e "\t* Pipeline found:\n\t\t${PIPE_ROOT}"

######################
# set up environment #
######################
echo "Setting up pipeline add-on environment."

# set project root directory
export PIPE_ADDON_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # allows startup.sh to be called from any directory.  will fail if last component of path to startup.sh is a symlink.

# set lmi reduction path
export REDUCTION="$PIPE_EXTRA_ROOT/lmi_reduction/"

# add pipeline directories to python path
export PYTHONPATH=$PIPE_EXTRA_ROOT/lmi_reduction:$PIPE_EXTRA_ROOT/lmi_observing:$PIPE_EXTRA_ROOT/lmi_photometry:$PYTHONPATH

echo "Pipeline add-on ready."