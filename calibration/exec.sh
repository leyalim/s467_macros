#!/bin/bash


function initialise () {
    #    for runnum in {237..380..1}
    for runnum in 238 273 361 366
    do
	list=$list' '$runnum
#	echo $runnum
    done
#    echo $list
}
    
#eval parallel --gnu --ungroup -j20 "root -l -b -q 'rawsofsci_offline.C('"{}"')'" ::: ${SEQ}
function myfunc () {
    #time nice root -l -b -q 'rawsofsci_offline.C('"$1"')' &> /dev/null
    time root -l -b -q 'tcal_VFTX_offline.C('"$1"')' &> /dev/null
    echo 'Finished run:'$1
}

function mapp() {
    #  This is from, http://prll.sourceforge.net/shell_parallel.html
    if [[ -z $MAPP_NR_CPUS ]] ; then
	#local MAPP_NR_CPUS=$(grep "processor:" < /proc/cpuinfo | wc -l)
	#   max core for calculation; modified by Toshiyuki Sumikama
	local MAPP_NR_CPUS=20
    fi
    local mapp_pid=$(exec bash -c 'echo $PPID')
    local mapp_funname=$1
    local -a mapp_params
    mapp_params=("$@")
    #   mapp_nr_args; modified by Toshiyuki Sumikama
    #  local mapp_nr_args=${#mapp_params[@]}
    local mapp_nr_args=`expr ${#mapp_params[@]} - 1`
    local mapp_current=0
    function mapp_trap() {
	#echo "MAPP PROGRESS: $((mapp_current*100/mapp_nr_args))%" 1>&2
	if [[ $mapp_current -lt $mapp_nr_args ]] ; then
	    let mapp_current+=1
	    (
		$mapp_funname "${mapp_params[$mapp_current]}"
		kill -USR1 $mapp_pid
	    ) &
	fi
    }

    trap mapp_trap SIGUSR1
    while [[ $mapp_current -lt $mapp_nr_args ]]; do
	wait
	if [[ $mapp_current -lt $mapp_nr_args && $? -lt 127 ]] ; then
	    sleep 1
	    local mapp_tmp_count=$mapp_current
	    wait
	    if [[ $mapp_tmp_count -eq $mapp_current ]] ; then
		echo "   MAPP_FORCE" 1>&2
		for i in $(seq 1 ${MAPP_NR_CPUS}) ; do
		    mapp_trap
		done
	    fi
	fi
    done
    for i in $(seq 1 ${MAPP_NR_CPUS}) ; do
	wait
    done
    trap - SIGUSR1
    unset -f mapp_trap
}


initialise
echo $list
(mapp myfunc $list)
unset list
echo 'done' #$list
