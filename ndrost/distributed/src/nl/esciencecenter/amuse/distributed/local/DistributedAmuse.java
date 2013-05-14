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
package nl.esciencecenter.amuse.distributed.local;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.octopus.Octopus;
import nl.esciencecenter.octopus.OctopusFactory;
import nl.esciencecenter.octopus.exceptions.OctopusException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Main Distributed AMUSE class. Started by AMUSE via the "Code" interface. Mostly contains objects that do the actual work.
 * 
 * @author Niels Drost
 * 
 */
public class DistributedAmuse {

    private static final Logger logger = LoggerFactory.getLogger(DistributedAmuse.class);

    //resources potentially available for starting reservations on
    private final ResourceManager resources;

    //starts pilots on resources. Also starts IPL server and hub when required
    private final ReservationManager reservations;

    //takes care of job queue, communicates with remote pilots
    private final JobManager jobs;

    //talks to AMUSE, handling any worker requests and messages
    private final WorkerConnectionServer workerConnectionHandler;

    //used to copy files, start jobs, etc.
    private final Octopus octopus;
    
    public DistributedAmuse() throws DistributedAmuseException {
        try {
            octopus = OctopusFactory.newOctopus(null);
        } catch (OctopusException e) {
            throw new DistributedAmuseException("could not create Octopus library object", e);
        }
        
        resources = new ResourceManager(octopus);

        reservations = new ReservationManager(octopus);

        jobs = new JobManager(reservations.getServerAddress());

        workerConnectionHandler = new WorkerConnectionServer(jobs);
    }

    public ResourceManager resources() {
        return resources;
    }

    public ReservationManager reservations() {
        return reservations;
    }

    public JobManager jobs() {
        return jobs;
    }

    /**
     * Port used by the IbisChannel to connect to when creating workers and stderr/stdout streams
     * 
     * @return the port used by the IbisChannel to connect to when creating workers and stderr/stdout streams
     */
    public int getWorkerPort() {
        logger.debug("returning worker port");
        return workerConnectionHandler.getPort();
    }

}
