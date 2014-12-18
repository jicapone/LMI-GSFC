#
#
#!/bin/bash

# check that RATIR pipeline has been sourced
[ -z "$RAT_PROJECT_ROOT" ] && echo "Need to run RATIR pipeline's startup.sh first.  Aborting." && exit 1

# set project root directory
export LMI_PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # allows startup.sh to be called from any directory.  will fail if last component of path to startup.sh is a symlink.

# set lmi reduction path
export REDUCTION="$LMI_PROJECT_ROOT/lmi_reduction/"

# add pipeline directories to python path
export PYTHONPATH=$LMI_PROJECT_ROOT/lmi_reduction:$LMI_PROJECT_ROOT/lmi_observing:$LMI_PROJECT_ROOT/lmi_photometry:$PYTHONPATH

# set pipeline aliases
alias cd_lmi='cd $LMI_PROJECT_ROOT'