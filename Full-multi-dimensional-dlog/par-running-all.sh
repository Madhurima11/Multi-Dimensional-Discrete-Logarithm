magma generate-target.m > out-generate-target.txt &
wait
magma pre-calculate-table.m > out-pre-calculate-table.txt &
magma multiplier-table-pre-calculate.m > out-multiplier-table-pre-calculate.txt &
wait
magma gs-full-log-finding.m > out-gs-full-log-finding.txt &
magma tt-full-log-finding.m > out-tt-full-log-finding.txt &
