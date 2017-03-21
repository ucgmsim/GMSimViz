try:
    from mpi4py import MPI
except ImportError:
    print "No mpi4py package available, exiting"
    exit(1)


ERROR = "ERROR"

class ParallelExecutor:
    def __init__(self):
        self.comm = MPI.COMM_WORLD
        self.cpu_number = self.comm.Get_size()

    def process_function_with_result(self, f, to_process, f_master = None):

        if len(to_process) == 0:
            return
        if self.cpu_number == 1:
            result = {0: []}
            # serial job
            for element in to_process:
                result[0].append(f(element))
            print "Terminating serial job"
            return result

        rank = self.comm.Get_rank()

        if rank == 0:
            currently_working = dict.fromkeys(range(self.cpu_number),"no")
            result = dict.fromkeys(range(self.cpu_number),[])

            worker_number = self.cpu_number
            initial_jobs = min(len(to_process)+1, worker_number)

            for i in range(1, initial_jobs):
                job = to_process.pop()
                self.comm.send(job, dest=i, tag=0)
                currently_working[i] = "yes"

            for i in range(initial_jobs, worker_number):
                self.comm.send(None, dest=i, tag=0)

            while True:
                status = MPI.Status()
                msg = self.comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)

                #free_worker, processor_result = self.comm.recv(source=MPI.ANY_SOURCE, tag=1)
                if status.Get_tag() == 1:
                    free_worker, processor_result = msg
                    # print processor_result
                    if processor_result == ERROR:
                        print "One task failed from processor %i" % free_worker
                    result[free_worker].append(processor_result)
                    currently_working[free_worker] = "no"
                    if len(to_process) > 0:
                        next_job = to_process.pop()
                        self.comm.send(next_job, dest=free_worker, tag=0)
                        currently_working[free_worker] = "yes"
                    else:
                        all_done = True
                        for processor in currently_working:
                            if currently_working[processor] == "yes":
                                all_done = False
                                continue

                        if all_done:
                            for i in range(1, initial_jobs):
                                # print "Sending shutdown to %i" % i
                                self.comm.send(None, dest=i, tag=0)
                            #print "Terminating master"
                            break
                #when tag == 2, let master execute f_master()
                elif status.Get_tag() == 2:
                    f_master(msg, status.Get_tag(), status.Get_source())

            return result

        else:
            job = self.comm.recv(source=0, tag=0)
            while not job is None:
                # process_function job
                try:
                    current_result = f(job)
                    # tell master that we are free
                    self.comm.send([rank, current_result], dest=0, tag=1)
                except Exception:
                    self.comm.send([rank, ERROR], dest=0, tag=1)

                job = self.comm.recv(source=0, tag=0)
                # print "%i received %s" % (rank,job)

            #print "Terminating worker %i" % rank
            #self.comm.send(rank, dest=0, tag=2)
            return None


def get_rank():
    return MPI.COMM_WORLD.Get_rank()


def get_barrier():
    return MPI.COMM_WORLD.Barrier()


def test_function(x):
    print "Dealing with %s" % str(x)


class Test:
    @staticmethod
    def test_method(x):
        print "Dealing with %s from class Test" % str(x)

    def non_static_method(self, x):
        return 1.0/x


if __name__ == "__main__":
    executor = ParallelExecutor()

    # easy_ints = range(1000)
    # executor.process_function(test_function, easy_ints)

    easy_ints = range(1000)
    t = Test()
    res = executor.process_function_with_result(t.non_static_method, easy_ints)
    if MPI.COMM_WORLD.Get_rank() == 0:
        values = res.values()
        result = []
        for value in values:
            result += value

        #print result
