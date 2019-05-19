ls|grep "x_"|awk '{print "cat "$1 "|grep \"~~~~~\" -A 22|awk '"'"'{print $2}'"'"'|tr '"'"'\\n'"'"' '"'"'\\t'"'"';echo"}'|bash
