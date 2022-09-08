#!/bin/bash

kerdir='../G_spool'

list_G_raw='list_G_raw'

ls -1 $kerdir | awk '{print "'$kerdir/'" $1}' > $list_G_raw
