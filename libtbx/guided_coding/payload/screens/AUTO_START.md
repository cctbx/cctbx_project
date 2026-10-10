AUTO START | <queue-id>
JOBS
<One line per job: id, what it fixes, and any dependency on another job.>
CHECKS
<The local checks each job must pass; the full server suite waits for the combined result.>
INSTALLATION
<The test-only installation used, and that working installations are not used.>
RESULTS
<Where saved candidates, tickets and the queue record will be. A usage limit stops the queue; resume is manual with /gc auto resume.>
STOP
<Say /gc stop here, or run this from any Terminal (the next line, exactly):>
python3 -I -B <absolute path of this procedure's auto_queue.py> stop
NOT AUTHORIZED
<Master and working installations are not changed; nothing is integrated or published; no new server access.>
ACTION: NONE NEEDED (reply stop to stop the queue)
