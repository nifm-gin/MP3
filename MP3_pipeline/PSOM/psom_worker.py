#!/usr/bin/env python
"""
This script starts a psom agent
"""
__author__ = 'poquirion'



import argparse
import os
import sys
import subprocess
import time


class Worker():
    def __init__(self, directory, worker_id):


        self.worker_id = worker_id
        self.process_done = False
        # only three types of worker
        if worker_id == 'manager':
            self.cmd = ['/bin/bash', '{0}/logs/tmp/psom_manager.sh'.format(directory, worker_id)]
            self.std_err = '{0}/logs/PIPE.err'.format(directory, worker_id)
            self.std_out = '{0}/logs/PIPE.out'.format(directory, worker_id)
            touch = ['PIPE.failed', 'PIPE.exit', 'PIPE.out']
            self.touch = ["{0}/logs/{1}".format(directory, t) for t in touch]
        elif worker_id == 'garbage':
            self.cmd = ['/bin/bash', '{0}/logs/tmp/psom_garbage.sh'.format(directory)]
            self.std_err = '{0}/logs/garbage/garbage.err'.format(directory)
            self.std_out = '{0}/logs/garbage/garbage.out'.format(directory)
            touch = ['garbage.failed', 'garbage.exit', 'garbage.out']
            self.touch = ["{0}/logs/garbage/{1}".format(directory, t) for t in touch]
        else:
            self.cmd =['/bin/bash', '{0}/logs/tmp/psom{1}.sh'.format(directory, worker_id)]
            self.std_err = '{0}/logs/worker/psom{1}/worker.err'.format(directory, worker_id)
            self.std_out = '{0}/logs/worker/psom{1}/worker.out'.format(directory, worker_id)
            touch = ['worker.failed', 'worker.exit', 'worker.out']
            self.touch = ["{0}/logs/worker/psom{1}/{2}".format(directory, worker_id, t) for t in touch]


        self.process = None
        self.fp_out = None
        self.fp_err = None

    def start(self):

        #Start agent
        self.fp_out = open(self.std_out, 'w')
        self.fp_err = open(self.std_err, 'w')

        print('execution {0}'.format(self.cmd))
        self.process = subprocess.Popen(self.cmd, stdout=self.fp_out, stderr=self.fp_err)

    def check_status(self):

        retcode = self.process.poll()
        # Let know manager that there was a problem
        if retcode is not None:
            print("psom{0} done with with returncode {1}".format(self.worker_id, retcode))
            for t in self.touch:
                print("touch {0}".format(t))
                with open(t, 'a'):
                    os.utime(t, None)
            self.fp_out.flush()
            self.fp_err.flush()
            self.process_done = True
        return retcode

def wait_till_all_process_are_done(workers, loop_sleep):
    """
    Return only when all workers have a process_done status to True

    :param workers: worker list
    :param loop_sleep: time to sleep after each loop
    :return:
    """

    all_done = False
    while not all_done:
        all_done = True
        for worker in workers:
            if not worker.process_done:
                worker.check_status()
                all_done = False
        time.sleep(loop_sleep)

    print("all workers are done with their jobs")

def main(args=None):

    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(description='Start a PSOM worker')

    parser.add_argument("--directory", "-d", type=str, required=True, help='The PSOM output directory')

    # parser.add_argument("--input_dir", "-i", type=str, required=False, help='The PSOM input directory')

    parser.add_argument("--worker_id", "-w", type=str, required=True, help='The PSOM given worker id')

    parsed = parser.parse_args(args)

    # We force the working directory to be at the root of the output directory
    os.chdir("{0}/..".format(parsed.directory))


    worker_per_node = os.getenv("PSOM_WORKER_PPN")
    if worker_per_node is None:
        worker_per_node = 1
    else:
        worker_per_node = int(worker_per_node)

    all_worker = []
    try:

        first = (int(parsed.worker_id) - 1) * worker_per_node + 1

        for worker_id in range(first, first+worker_per_node):
            w = Worker(parsed.directory, worker_id)
            w.start()
            all_worker.append(w)

    except ValueError:
        w = Worker(parsed.directory, parsed.worker_id)
        w.start()
        all_worker.append(w)


    wait_till_all_process_are_done(all_worker, 5)

if __name__ == '__main__':
    # main(["-d", "/home/poquirion/simexp/test/result", "-w", "manager"])
    main()
