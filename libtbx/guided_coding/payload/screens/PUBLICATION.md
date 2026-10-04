PUBLICATION | <batch-id> | <manifest-sha256>
BATCH
<Exactly what would be pushed and its target: each outgoing commit from OUTGOING.txt in full, its branch and remote; nothing else moves and no tag is created.>
SERVER CHECK
<Full-suite result and roster compared with baseline, referencing packet logs; or SERVER_SUITE: NOT RUN with the Developer's quoted waiver for this exact batch.>
LIMITS
<Material gaps, any accepted risk, whether the exact combined tree had a full run, and which installations the plan names for a separate update (publication authorizes none).>
GATE
Reading SHA256: <sha256-of-reading-file>
Verdict: <PROCEED or PROCEED IF condition; if conditional add Disposition SHA256: <sha256> after this line>
DECISION
<Publish pushes this batch after the precheck, with the explicit refspec and tag following disabled; hold keeps it local. The reading above is scoped to publication or to both gates of this packet.>
ACTION: PUBLISH / HOLD
