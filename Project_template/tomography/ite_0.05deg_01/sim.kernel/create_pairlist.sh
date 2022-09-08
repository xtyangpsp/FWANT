#!/bin/bash

cat *_conf | grep "/" | awk '{if(NF==1) print $1}'| awk -F/ '{print $(NF-3), $(NF-2)}'| sort | uniq > .temp


