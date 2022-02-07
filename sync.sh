#!/bin/sh

#  sync.sh
#  
#
#  Created by Andrzej Chlebicki on 23.12.2018.
#  


# Sources
rsync -r -avz --prune-empty-dirs --exclude '*.html' --exclude '*.csv' --exclude '*.png' --exclude '.*' --exclude "tasks/*" --exclude 'Z*_*/*' --exclude 'ON*_*/*' --exclude 'lowd*_*/*' --exclude 'tricritical*_*/*' --exclude 'ising*_*/*' --exclude 'notatki' --exclude 'archive' . ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalRenormalization

# Results
rsync -r -avz --exclude '*.sh*' --exclude "tasks/*" --exclude '*.wl' --exclude '*.nb' --exclude 'flow_equations' ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalRenormalization/ .
