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
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.local.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.local.Reservation;
import nl.esciencecenter.amuse.distributed.local.Resource;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Interface from generated code interface to the distributed amuse implementation. Simply forwards all calls. Must be in the
 * "default" package as the CodeInterface and Worker are also generated there.
 * 
 * @author Niels Drost
 * 
 */
public class Code implements CodeInterface {

    private static final Logger logger = LoggerFactory.getLogger(Code.class);

    private final DistributedAmuse distributedAmuse;

    public Code() throws DistributedAmuseException {
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
    public int new_resource(int[] index_of_the_resource, String[] name, String[] hostname, String[] amuse_dir, int[] port,
            String[] username, String[] scheduler_type, int[] start_hub, int count) {
        try {
            for (int i = 0; i < count; i++) {
                Boolean startHub;
                
                if (start_hub[i] == -1) {
                    //default: auto
                    startHub = null;
                } else if (start_hub[i] == 0) {
                    startHub = false;
                } else {
                    startHub = true;
                }
                Resource resource =
                        distributedAmuse.resources().newResource(name[i], hostname[i], amuse_dir[i], port[i], username[i], scheduler_type[i], startHub);
                index_of_the_resource[i] = resource.getId();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on running distributed code: " + e);
            return 10;
        }
    }

    @Override
    public int delete_resource(int[] index_of_the_resource, int count) {
        try {
            for (int i = 0; i < count; i++) {
                Resource resource = distributedAmuse.resources().getResource(index_of_the_resource[i]);
                distributedAmuse.resources().deleteResource(resource);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on running distributed code: " + e);
            return 10;
        }
    }

    @Override
    public int new_reservation(int[] reservation_id, String[] resource_name, String[] queue_name, int[] node_count,
            int[] time_minutes, String[] node_label, int count) {
        try {
            for (int i = 0; i < count; i++) {
                Resource resource = distributedAmuse.resources().getResource(resource_name[i]);
                Reservation result =
                        distributedAmuse.reservations().newReservation(resource, queue_name[i], node_count[i],
                                time_minutes[i], node_label[i]);

                reservation_id[i] = result.getID();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on running distributed code: " + e);
            return 10;
        }

    }

    @Override
    public int delete_reservation(int[] reservation_id, int count) {
        try {
            for (int i = 0; i < count; i++) {
                distributedAmuse.reservations().deleteReservation(reservation_id[i]);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on running distributed code: " + e);
            return 10;
        }
    }

    @Override
    public int wait_for_reservations() {
        distributedAmuse.reservations().waitForAllReservations();
        return 0;
    }

    @Override
    public int submit_pickled_function_job(int[] job_id, String[] function, String[] arguments, String[] node_label, int count) {
        try {
            for (int i = 0; i < count; i++) {
                job_id[i] = distributedAmuse.jobs().submitPickledJob(function[i], arguments[i], node_label[i]);

            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on running distributed code: " + e);
            return 10;
        }

    }

    @Override
    public int submit_script_job(int[] job_id, String[] script, String[] arguments, String[] code_dir, String[] node_label,
            int[] re_use_code_files, int count) {
        try {
            for (int i = 0; i < count; i++) {
                boolean useCodeCache = re_use_code_files[i] != 0;
                job_id[i] =
                        distributedAmuse.jobs().submitScriptJob(script[i], arguments[i], code_dir[i], node_label[i],
                                useCodeCache);

            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on running distributed code: " + e);
            return 10;
        }

    }

    @Override
    public int get_pickled_function_job_result(int[] job_id, String[] result, int count) {
        try {
            for (int i = 0; i < count; i++) {
                result[i] = distributedAmuse.jobs().getJobResult(job_id[i]);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on running distributed code: " + e);
            return 10;
        }
    }

    @Override
    public int wait_for_jobs() {
        try {
            distributedAmuse.jobs().waitForAllJobs();
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on running distributed code: " + e);
            return 10;
        }

    }

}
