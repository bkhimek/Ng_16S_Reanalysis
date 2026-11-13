#!/usr/bin/env bash
# status.sh - is a pipeline step running, finished or failed?
# Usage (repository root):  bash scripts/status.sh [logfile]
# Without an argument it checks the most recently modified *.out log under
# results/ and data/raw/.
set -uo pipefail

log="${1:-$(ls -t results/*.out results/*/*.out data/raw/*.out 2>/dev/null | head -1)}"
[[ -n "$log" && -f "$log" ]] || { echo "No log file found."; exit 1; }

procs=$(pgrep -af "scripts/[0-9][0-9]_|qiime |R --(no-echo|slave)|Rscript|run_dada" | grep -v "status.sh" || true)
age=$(( $(date +%s) - $(stat -c %Y "$log") ))
last_step=$(grep -E "^(\+ |=== |DADA2 arm)" "$log" | tail -1 | cut -c1-140)
errors=$(grep -nE "Traceback|^Error|Plugin error|There were some problems|No such option|FAILED|Killed|MemoryError|No space left" "$log" | tail -3)

if [[ -n "$procs" ]]; then
  verdict="RUNNING"
elif [[ -n "$errors" ]]; then
  verdict="FAILED"
elif tail -3 "$log" | grep -qE "^Done|verified\.$"; then
  verdict="FINISHED"
else
  verdict="STOPPED (no process running and no 'Done' line: check the log)"
fi

echo "Log:            $log"
echo "Status:         $verdict"
echo "Current step:   ${last_step:-n/a}"
printf "Log last changed: %dh %02dm ago\n" $((age/3600)) $(((age%3600)/60))
if [[ -n "$procs" ]]; then
  echo "Processes:"
  echo "$procs" | cut -c1-110 | sed 's/^/  /'
  echo "CPU of the busiest process:"
  ps -eo pcpu,etime,comm --sort=-pcpu | sed -n 2p | sed 's/^/  %CPU, elapsed, name: /'
fi
if [[ -n "$errors" ]]; then
  echo "Error lines:"; echo "$errors" | cut -c1-140 | sed 's/^/  /'
fi
echo "Last 3 log lines:"; tail -3 "$log" | cut -c1-140 | sed 's/^/  /'
