#!/bin/sh

RNAMOTIFSCANX_PATH=$1
echo "Setting RNAMOTIFSCANX_PATH..."
export RNAMOTIFSCANX_PATH
echo "Setting RNAVIEW_PATH..."
RNAVIEW=$RNAMOTIFSCANX_PATH/thirdparty/RNAVIEW 
export RNAVIEW
#echo "Trying to compile RNAVIEW..."
#cd $RNAMOTIFSCANX_PATH/thirdparty/RNAVIEW
#make
#cd $RNAMOTIFSCANX_PATH
echo "Setting environment done."
