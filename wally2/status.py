#!/usr/bin/env python3
import sys
import subprocess
import logging

# Max number of attempts
MAX_ATTEMPTS = 100

# Define string constants for statuses
failed = 'failed'
success = 'success'
running = 'running'

# Setup logger (check if this is needed)
logger = logging.getLogger('__name__')

# Dictionary giving the snakemake state for a job for each SLURM status
status_table = {'BOOT_FAIL': failed,
                'CANCELLED': failed,
                'COMPLETED': success,
                'CONFIGURING': running,
                'COMPLETING': running,
                'DEADLINE': failed,
                'FAILED': failed,
                'NODE_FAIL': failed,
                'OUT_OF_MEMORY': failed,
                'PENDING': running,
                'PREEMPTED': failed,
                'RUNNING': running,
                'RESV_DEL_HOLD': running,
                'REQUEUE_FED': failed,
                'REQUEUE_HOLD': failed,
                'REQUEUED': failed,
                'RESIZING': failed,
                'REVOKED': failed,
                'SIGNALING': failed,
                'SPECIAL_EXIT': failed,
                'STAGE_OUT': failed,
                'STOPPED': failed,
                'SUSPENDED': failed,
                'TIMEOUT': failed}

# Recover job ID from CL arguments
job_id = int(sys.argv[1])

# Generate status checking command from job ID
cmd = 'sacct -nbPj {}'.format(job_id)

# Repeat status check until successful or until MAX_ATTEMPTs is reached
status = None
for attempt in range(MAX_ATTEMPTS):
    try:
        status_string = subprocess.check_output(cmd, shell=True)  # Run command to get job status
        status = status_string.decode('utf-8').split('|')[1]  # Response is "ID|STATUS|EXIT CODE"
        break  # Command was successful
    except subprocess.CalledProcessError as error:
        logger.error('Failed status checking attempt {attempt}/{MAX_ATTEMPTS}: <{error}>'.format(attemt=attempt, MAX_ATTEMPTS=MAX_ATTEMPTS, error=error))
    except IndexError as error:
        # Split failed, meaning the response was not following the expected format
        # This happens for instance if the job id does not exist
        logger.error('Unexpected response in checking attempt {attempt}/{MAX_ATTEMPTS}: <{status}>'.format(attemt=attempt, MAX_ATTEMPTS=MAX_ATTEMPTS, status=status_string.decode("utf-8")))
        pass

# Prints the corresponding snakemake status
if status not in status_table:
    logger.error('Unknown status: <{}>'.format(status))
    print('failed')  # If something unexpected happened, the job is marked as failed by default
else:
    print(status_table[status])
