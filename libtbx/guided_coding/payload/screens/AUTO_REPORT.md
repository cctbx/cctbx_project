AUTO REPORT | <queue-id>
JOBS
<The table printed by auto_queue.py report, copied exactly, with no added or changed rows.>
IN SIMPLE WORDS
<What is ready, what is blocked and why, what is still waiting, in ordinary words.>
Tickets: <usable path to each ready job's ticket>
APPROVAL SUMMARIES
<For each ready job, in table order, seven contiguous lines (no blank lines between blocks): 'Job' and the job's own id with a colon, then its title (for example 'Job A: launcher fix'), then 'Bug: ', 'Fix: ', 'Test: ' (failed before, passes after), 'Criterion: ' (chosen by the Guide), 'Limits: ' (material limits and any behaviour change), 'Approval: accepts this exact packet and merges nothing.' Or the single line 'None ready.'>
NEXT
<How to answer in one message, for example: approve A C; revise B: what to change; discard D.>
Recommendation: <which jobs to approve or inspect, and why.>
ACTION: APPROVE / INSPECT / REVISE / DISCARD
