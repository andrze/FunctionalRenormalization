#!/bin/sh

#  sync.sh
#  
#
#  Created by Andrzej Chlebicki on 23.12.2018.
#  


# Sources
rsync -r -avz --prune-empty-dirs --exclude '*.html' --exclude '*.csv' --exclude '*.png' --exclude '.*' --exclude 'Z*_*/*' --exclude 'ON*_*/*' . ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalRenormalization

# Results
rsync -r -avz --exclude '*.sh*' --exclude '*.wl' --exclude '*.nb' --exclude 'flow_equations' ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalRenormalization/ .
