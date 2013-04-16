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
package nl.esciencecenter.amuse.distributed;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DistributedAmuse {

    private static final Logger logger = LoggerFactory.getLogger(DistributedAmuse.class);

    private final int workerPort = 6644;
    
    /**
     * Port used by the IbisChannel to connect to when creating workers and stderr/stdout streams
     * 
     * @return the port used by the IbisChannel to connect to when creating workers and stderr/stdout streams
     */
    public int getWorkerPort() {
        logger.debug("returning worker port " + workerPort);
        return workerPort;
    }

    public int newResource(String name, String hostname, String amuseDir, int port, String username, String schedulerType) {
        logger.debug("creating new resource: name = " + name + " hostname = " + hostname + " port = " + port + " user name = "
                + username + " scheduler type = " + schedulerType + " amuse dir = " + amuseDir);
        return 0;
    }

    public void deleteResource(int resourceID) {
        logger.debug("deleting resource " + resourceID);

        // TODO Auto-generated method stub

    }

    public int newReservation(String resourceName, String queueName, int nodeCount, int timeMinutes, String nodeLabel) {
        logger.debug("reserving new nodes: resource name = " + resourceName + " queue name = " + queueName
                + " number of nodes = " + nodeCount + " time (in minutes) = " + timeMinutes + " node label = " + nodeLabel);

        return 0;
    }

    public void deleteReservation(int reservationID) {
        logger.debug("deleting reservation " + reservationID);

    }

    public int submitJob(PickledJobDescription jobDescription) {
        logger.debug("submitting new job: " + jobDescription);

        return 0;
    }

    public int submitJob(ScriptJobDescription jobDescription) {
        logger.debug("submitting new job: " + jobDescription);
        return 0;
    }

    public String getJobResult(int jobID) {
        logger.debug("gettig job result for " + jobID);

        // TODO Auto-generated method stub
        return "hurray job done!";
    }

    public void waitForAllJobs() {
        logger.debug("waiting for all jobs to complete");

        // TODO Auto-generated method stub

    }

    public void waitForAllReservations() {
        logger.debug("waiting for all reservations to start");

        // TODO Auto-generated method stub

    }

}
