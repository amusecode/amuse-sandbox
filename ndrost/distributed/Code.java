/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.PickledJobDescription;
import nl.esciencecenter.amuse.distributed.ScriptJobDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Code implements CodeInterface {
    
    private final DistributedAmuse distributedAmuse;

    public Code() {
        distributedAmuse = new DistributedAmuse();
    }

    @Override
    public int initialize_code() {
        //IGNORED
        return 0;
    }
    
    @Override
    public int commit_parameters() {
        //IGNORED
        return 0;
    }
    
    @Override
    public int recommit_parameters() {
        //IGNORED
        return 0;
    }
    
    @Override
    public int get_port() {
        return distributedAmuse.getWorkerPort();
    }

    @Override
    public int new_resource(int[] index_of_the_resource, String[] name, String[] hostname, String[] amuse_dir, int[] port, String[] username, String[] scheduler_type, int count) {
        for (int i = 0; i < count; i++) {
            index_of_the_resource[i] =
                    distributedAmuse.newResource(name[i], hostname[i],  amuse_dir[i], port[i], username[i], scheduler_type[i]);
        }
        return 0;
    }

    @Override
    public int delete_resource(int[] index_of_the_resource, int count) {
        for (int i = 0; i < count; i++) {
            distributedAmuse.deleteResource(index_of_the_resource[i]);
        }
        return 0;
    }

    @Override
    public int new_reservation(int[] reservation_id, String[] resource_name, String[] queue_name, int[] node_count,
            int[] time_minutes, String[] node_label, int count) {
        for (int i = 0; i < count; i++) {
            reservation_id[i] =
                    distributedAmuse.newReservation(resource_name[i], queue_name[i], node_count[i], time_minutes[i],
                            node_label[i]);
        }
        return 0;
    }

    @Override
    public int delete_reservation(int[] reservation_id, int count) {
        for (int i = 0; i < count; i++) {
            distributedAmuse.deleteReservation(reservation_id[i]);
        }
        return 0;
    }

    @Override
    public int submit_pickled_function_job(int[] job_id, String[] function, String[] arguments, String[] node_label, int count) {
        for(int i = 0; i < count; i++) {
            job_id[i] = distributedAmuse.submitJob(new PickledJobDescription(function[i], arguments[i], node_label[i]));
        }
        return 0;
    }

    @Override
    public int submit_script_job(int[] job_id, String[] script, String[] arguments, String[] code_dir, String[] node_label,
            int[] re_use_code_files, int count) {
        for(int i = 0; i < count; i++) {
            boolean useCodeCache = re_use_code_files[i] != 0;
            job_id[i] = distributedAmuse.submitJob(new ScriptJobDescription(script[i], arguments[i], code_dir[i], node_label[i], useCodeCache));
        }
        return 0;
    }

    @Override
    public int get_pickled_function_job_result(int[] job_id, String[] result, int count) {
        // TODO Auto-generated method stub
        for (int i = 0; i < count; i++) {
            result[i] = distributedAmuse.getJobResult(job_id[i]);
        }

        return 0;
    }

    @Override
    public int wait_for_jobs() {
          distributedAmuse.waitForAllJobs();
          return 0;
    }


    @Override
    public int wait_for_reservations() {
        distributedAmuse.waitForAllReservations();
        return 0;
    }
}
