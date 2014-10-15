#
#
#!/bin/bash

PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # allows startup.sh to be called from any directory.  will fail if last component of path to startup.sh is a symlink.

# add pipeline directories to python path
export PYTHONPATH=$PYTHONPATH:$PROJECT_ROOT/lmi_reduction:$PROJECT_ROOT/lmi_observing:$PROJECT_ROOT/lmi_photometry
