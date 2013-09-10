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
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.reservations.Reservation;
import nl.esciencecenter.amuse.distributed.resources.Resource;

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

    public Code(String codeDir, String amuseRootDir) throws DistributedAmuseException {
        distributedAmuse = new DistributedAmuse(codeDir, amuseRootDir, 0);
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
    public int get_worker_port() {
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
                        distributedAmuse.resourceManager().newResource(name[i], hostname[i], amuse_dir[i], port[i], username[i],
                                scheduler_type[i], startHub);
                index_of_the_resource[i] = resource.getId();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on creating new resource: " + e, e);
            return 10;
        }
    }

    @Override
    public int delete_resource(int[] index_of_the_resource, int count) {
        try {
            for (int i = 0; i < count; i++) {
                Resource resource = distributedAmuse.resourceManager().getResource(index_of_the_resource[i]);
                distributedAmuse.resourceManager().deleteResource(resource);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on deleting resource: " + e, e);
            return 10;
        }
    }

    @Override
    public int new_reservation(int[] reservation_id, String[] resource_name, String[] queue_name, int[] node_count,
            int[] time_minutes, String[] node_label, int count) {
        try {
            for (int i = 0; i < count; i++) {
                Reservation result =
                        distributedAmuse.reservationManager().newReservation(resource_name[i], queue_name[i], node_count[i],
                                time_minutes[i], node_label[i]);

                reservation_id[i] = result.getID();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on creating new reservation: " + e, e);
            return 10;
        }

    }

    @Override
    public int delete_reservation(int[] reservation_id, int count) {
        try {
            for (int i = 0; i < count; i++) {
                distributedAmuse.reservationManager().deleteReservation(reservation_id[i]);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on deleting reservation: " + e, e);
            return 10;
        }
    }

    @Override
    public int wait_for_reservations() {
        try {
            distributedAmuse.reservationManager().waitForAllReservations();
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on waiting for reservations: " + e, e);
            return 10;
        }
    }

    @Override
    public int submit_pickled_function_job(int[] job_id, String[] function, String[] arguments, String[] node_label, int count) {
        try {
            for (int i = 0; i < count; i++) {
                job_id[i] = distributedAmuse.jobManager().submitPickledJob(function[i], arguments[i], node_label[i]);

            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Cannot submit pickled job: " + e, e);
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
                        distributedAmuse.jobManager().submitScriptJob(script[i], arguments[i], code_dir[i], node_label[i],
                                useCodeCache);

            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Cannot submit script job: " + e, e);
            return 10;
        }

    }

    @Override
    public int get_pickled_function_job_result(int[] job_id, String[] result, int count) {
        try {
            for (int i = 0; i < count; i++) {
                result[i] = distributedAmuse.jobManager().getJobResult(job_id[i]);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on getting job result: " + e, e);
            return 10;
        }
    }

    @Override
    public int wait_for_jobs() {
        try {
            distributedAmuse.jobManager().waitForAllBatchJobs();
            return 0;
        } catch (DistributedAmuseException e) {
            logger.error("Error on waiting for jobs: " + e, e);
            return 10;
        }

    }

    /**
     * 
     */
    @Override
    public void end() {
        distributedAmuse.end();
    }

}
