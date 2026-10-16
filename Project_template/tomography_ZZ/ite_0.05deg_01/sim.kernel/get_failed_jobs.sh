sacct | grep FAILED | grep bell | grep kn | awk '{print $2}' | sed -e 's/kn/submit/g'
