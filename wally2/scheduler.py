import logging
import os
import re
import subprocess
import sys
import time
import yaml
from collections import defaultdict
from snakemake.utils import read_job_properties


def convert_time(runtime):
    '''
    Utility function to convert a runtime from SLURM runtime format to seconds.
    Acceptable time formats for SLURM:
    - "minutes"
    - "minutes:seconds"
    - "hours:minutes:seconds"
    - "days-hours"
    - "days-hours:minutes"
    - "days-hours:minutes:seconds"
    '''
    runtime = str(runtime)
    d, h, m, s = 0, 0, 0, 0
    fields = runtime.split('-')
    parsing_error = False
    if len(fields) == 2:
        d = int(fields[0])
        fields = tuple(int(f) for f in fields[1].split(':'))
        if len(fields) == 1:
            h = fields[0]
        elif len(fields) == 2:
            h, m = fields
        elif len(fields) == 3:
            h, m, s = fields
        else:
            parsing_error = True
    elif len(fields) == 1:
        fields = tuple(int(f) for f in fields[0].split(':'))
        if len(fields) == 1:
            m = fields[0]
        elif len(fields) == 2:
            m, s = fields
        elif len(fields) == 3:
            h, m, s = fields
        else:
            parsing_error = True
    else:
        parsing_error = True
    if parsing_error:
        #logging.error(f'Invalid format for runtime string <{runtime}> (accepted formats: "M", "M:S", "H:M:S", "D-H", "D-H:M", "D-H:M:S"')
        logging.error('Invalid format for runtime string <{}> (accepted formats: "M", "M:S", "H:M:S", "D-H", "D-H:M", "D-H:M:S"'.format(runtime))
        exit(1)
    else:
        return ((((d * 24) + h) * 60) + m) * 60 + s


def convert_time_slurm(runtime):
    '''
    Utility function to convert a runtime from seconds to format 'D-HH:MM:SS'
    '''
    d, h, m, s = (0, 0, 0, 0)
    s = runtime % 60
    runtime = (runtime - s) // 60
    m = runtime % 60
    runtime = (runtime - m) // 60
    h = runtime % 24
    runtime = (runtime - h) // 24
    d = runtime
    if h < 10:
        #h = f'0{h}'
        h = '0{}'.format(h)
    if m < 10:
        #m = f'0{m}'
        m = '0{}'.format(m)
    if s < 10:
        #s = f'0{s}'
        s = '0{}'.format(s)
    #return f'{d}-{h}:{m}:{s}'
    return '{}-{}:{}:{}'.format(d,h,m,s)


def output(cmd):
    '''
    Wrapper around subprocess.check_output that returns the output in utf-8 format
    '''
    return subprocess.check_output(cmd, shell=True).decode('utf-8')


class SlurmScheduler:

    def __init__(self):
        self.cfg = {}
        self.partitions_info = {}
        self.submission_settings = defaultdict(lambda: None)
        self.command = ''
        self.partitions_file = ''
        self.jobscript = sys.argv[1]
        self.job_properties = read_job_properties(self.jobscript)
        self.load_slurm_config()
        self.update_partitions_info()
        self.load_partitions_info()

    def load_slurm_config(self, cfg_path='slurm.yaml'):
        '''
        Loads settings from the slurm config YAML file into a dictionary.
        '''
        cfg_path = os.path.join(os.path.split(os.path.abspath(__file__))[0], cfg_path)
        self.cfg = yaml.safe_load(open(cfg_path))
        if not self.cfg['blacklist']:
            self.cfg['blacklist'] = set()
        self.partitions_file = os.path.join(os.path.split(os.path.abspath(__file__))[0], self.cfg['scheduler']['partitions_file'])

    def update_partitions_info(self):
        '''
        Retrieve information about partitions available on the cluster and store the data in a yaml file.
        This information is updated at a rate specified in the 'slurm.yaml' config file.
        Only partitions on which the user is allowed to submit are retained.
        '''
        if os.path.isfile(self.partitions_file):
            time_since_modification = time.time() - os.path.getmtime(self.partitions_file)
            if time_since_modification < self.cfg['scheduler']['partitions_update_days'] * 24 * 60 * 60:
                return  # Only update file if it's not there or it's older than update rate
        username = output('whoami').strip()  # Retrieve username
        # Retrieve slurm account
        # The command's response has format "<user>|<account>|<Admin>".
        # The account name is retrieved from the second field
        #account = output(f'sacctmgr -Pn show user {username}').split('|')[1].strip()
        account = output('sacctmgr -Pn show user {}'.format(username)).split('|')[1].strip()
        # Retrieve groups
        # The command's response has format "<user> : <group1> <group2> ... <groupN>"
        # The response is parsed to create a list of groups
        #groups = [g for g in output(f'groups {username}')[:-1].split(':')[-1].split(' ') if g != '']
        #groups = [g for g in output('groups {}'.format(username))[:-1].split(':')[-1].split(' ') if g != '']
        groups = [g for g in curGroups[:-1].split(':')[-1].split(' ') if g != '']
        # Retrieve and parse partition information
        # The command's response is a variable space-separated table with a header and a line per partition's
        # node configuration (CPU, memory ...). This means there are usually multiple lines for each unique partition.
        # The response is parsed to create a list in which each element is a list of values for a partition's node configuration,
        # with the first element being field names given by the header
        sinfo_response = output('sinfo --noconvert -eO "partitionname,cpus,memory,time,maxcpuspernode,groups,available,prioritytier"')
        info = [[field for field in partition[:-1].split(' ') if field != ''] for
                partition in sinfo_response.split('\n') if partition != '']
        # Converts the list of list into a list of dictionaries, with field names given by the header as keys and settings for
        # the given partition as values.
        info = [{k: v for k, v in zip(info[0], partition)} for partition in info[1:]]
        # Final partition info format: dict with partition name as key and dict of settings as values
        summary = defaultdict(lambda: defaultdict(lambda: None))
        # Parse output from sinfo to create the final partition information summary. When there are multiple node configurations
        # for a partition, the maximum values are retained
        for partition in info:
            name = partition['PARTITION']
            # Use scontrol to get list of accounts allowed to submit to current partition
            #scontrol_response = output(f'scontrol show partition {name} | grep AllowAccounts')
            scontrol_response = output('scontrol show partition {} | grep AllowAccounts'.format(name))
            allowed_accounts = re.search(r'AllowAccounts=([^\s]+)', scontrol_response).group(1)
            # Check that current user's account is allowed to submit, otherwise discards partition
            if allowed_accounts != 'ALL' and account not in allowed_accounts.split(','):
                continue
            for field, value in partition.items():
                if value.isdigit():  # For integer values, replace if current value is bigger than saved value
                    if not summary[name][field] or int(value) > summary[name][field]:
                        summary[name][field] = int(value)
                elif field == 'TIMELIMIT' and value != 'infinite':
                    duration = convert_time(value)  # Convert max runtime to seconds for easy comparison
                    if not summary[name][field] or duration > summary[name][field]:  # Replace if bigger
                        summary[name][field] = duration
                else:
                    summary[name][field] = value
        # Export summary to a yaml file
        with open(self.partitions_file, 'w') as summary_file:
            yaml.dump(dict({k: dict(v) for k, v in summary.items()}), summary_file, default_flow_style=False)

    def load_partitions_info(self):
        '''
        Load partition info YAML file generated by 'update_partitions_info' into a dictionary
        Filter out partitions blacklisted in 'slurm.yaml' config file
        '''
        self.partitions_info = yaml.safe_load(open(self.partitions_file))
        self.partitions_info = {k: v for k, v in self.partitions_info.items() if k not in self.cfg['blacklist']}

    def check_for_setting(self, setting):
        '''
        Check if a given setting is present in the job's properties (top-level setting)
        '''
        return self.job_properties[setting] if setting in self.job_properties else None

    def check_for_param(self, param):
        '''
        Check if a given setting is present in the job's parameters
        '''
        return self.job_properties['params'][param] if 'params' in self.job_properties and param in self.job_properties['params'] else None

    def check_for_resource(self, resource):
        '''
        Check if a given setting is present in the job's resources
        '''
        return self.job_properties['resources'][resource] if 'resources' in self.job_properties and resource in self.job_properties['resources'] else None

    def get_submission_settings(self):
        '''
        For each submission setting defined in the 'slurm.yaml' <options> field, check if a value is specified
        in the job's properties. First check if setting is defined in 'resources', then in 'params', and then in
        top-level settings.
        '''
        for setting, arg_string in self.cfg['options'].items():
            value = self.check_for_resource(setting)
            if value is None:
                value = self.check_for_param(setting)
            if value is None:
                value = self.check_for_setting(setting)
            if value is not None and value != [] and value != {}:
                if setting == 'log':
                    value = value[0]
                    dir_path = os.path.abspath(os.path.split(value)[0])
                    if not os.path.isdir(dir_path):
                        os.makedirs(dir_path)
                self.submission_settings[setting] = value
        if 'runtime' in self.submission_settings:  # Runtime was specified in slurm format, convert to seconds for comparison
            self.submission_settings['runtime'] = convert_time(self.submission_settings['runtime'])
        elif 'runtime_s' in self.submission_settings:  # Runtime was specified in seconds, update runtime value with runtime_s
            self.submission_settings['runtime_s'] = int(self.submission_settings['runtime_s'])

    def set_partition(self):
        '''
        If a partition was not manually specified by the user in the snakefile, look for the best partition
        according to the job's resources requirements. At the moment, partitions are ranked based only on PRIO_TIER,
        prioritizing partitions with higher tiers. For each partition in order of desirability, the function checks
        if the partition is up and if threads, memory, and runtime requirements can be satisfied. The first partition
        to satisfy these criteria is used. If no suitable partition is found, the function prints an error and exits.
        '''
        if self.submission_settings['partition'] is None:
            for partition, data in sorted(self.partitions_info.items(), key=lambda x: x[1]['PRIO_TIER'], reverse=True):
                suitable = True
                if data['AVAIL'] != 'up':
                    suitable = False
                elif self.submission_settings['threads'] and int(self.submission_settings['threads']) > int(data['CPUS']):
                    suitable = False
                elif self.submission_settings['memory'] and int(self.submission_settings['memory']) > int(data['MEMORY']):
                    suitable = False
                elif self.submission_settings['mem_mb'] and int(self.submission_settings['mem_mb']) > int(data['MEMORY']):
                    suitable = False
                elif self.submission_settings['runtime'] and self.submission_settings['runtime'] > int(data['TIMELIMIT']):
                    suitable = False
                elif self.submission_settings['runtime_s'] and self.submission_settings['runtime_s'] > int(data['TIMELIMIT']):
                    suitable = False
                if suitable:
                    self.submission_settings['partition'] = partition
                    return
            logging.error('No partition was found to satisfy resources requirements')
            exit(1)
        elif self.submission_settings['partition'] not in self.partitions_info:
            #logging.error(f'Partition <{self.submission_settings["partition"]}> specified by user was not found')
            logging.error('Partition <{}> specified by user was not found'.format(self.submission_settings["partition"]))
            exit(1)

    def generate_command(self):
        '''
        Generate the submission command based on submission settings.
        '''
        if self.submission_settings['runtime']:
            self.submission_settings['runtime'] = convert_time_slurm(self.submission_settings['runtime'])
        if self.submission_settings['runtime_s']:
            self.submission_settings['runtime_s'] = convert_time_slurm(self.submission_settings['runtime_s'])
        self.command = 'sbatch '
        for setting, arg_string in self.cfg['options'].items():
            if self.submission_settings[setting] is not None:
                self.command += arg_string.format(*([str(self.submission_settings[setting])] * arg_string.count('{}'))) + ' '
        #self.command += f'{self.jobscript}'
        self.command += '{}'.format(self.jobscript)

    def submit_command(self):
        '''
        Submit the command and parse response to output only the submission id, which is required by Snakemake to check
        the job's status.
        '''
        submit_response = subprocess.run(self.command, check=True, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
        submit_regex = re.search(r'Submitted batch job (\d+)', submit_response)
        print(submit_regex.group(1))

    def submit(self):
        '''
        Top-level function handling all steps of the submission process.
        '''
        self.get_submission_settings()
        self.set_partition()
        self.generate_command()
        self.submit_command()
