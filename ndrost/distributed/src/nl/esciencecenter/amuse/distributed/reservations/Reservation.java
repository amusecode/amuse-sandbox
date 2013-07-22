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
package nl.esciencecenter.amuse.distributed.reservations;

import java.net.URI;
import java.util.List;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.resources.Resource;
import nl.esciencecenter.octopus.Octopus;
import nl.esciencecenter.octopus.credentials.Credential;
import nl.esciencecenter.octopus.exceptions.OctopusException;
import nl.esciencecenter.octopus.exceptions.OctopusIOException;
import nl.esciencecenter.octopus.jobs.Job;
import nl.esciencecenter.octopus.jobs.JobDescription;
import nl.esciencecenter.octopus.jobs.Scheduler;
import nl.esciencecenter.octopus.util.JavaJobDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class Reservation {

    private static final Logger logger = LoggerFactory.getLogger(Reservation.class);

    private static int nextID = 0;

    private static int getNextID() {
        return nextID++;
    }

    private static JavaJobDescription createJobDesciption(Resource resource, String queueName, int nodeCount, int timeMinutes,
            String nodeLabel, String serverAddress, String[] hubAddresses) throws DistributedAmuseException {
        JavaJobDescription result = new JavaJobDescription();

        result.setInteractive(false);
        
        if (queueName != null && !queueName.isEmpty()) {
            result.setQueueName(queueName);
        }
        
        result.setNodeCount(nodeCount);
        result.setMaxTime(timeMinutes);

        AmuseConfiguration configuration = resource.getConfiguration();

        result.setExecutable(configuration.getJava());

        List<String> classpath = result.getJavaClasspath();
        classpath.add(configuration.getAmuseHome().getPath() + "/sandbox/ndrost/distributed/distributed-server.jar");
        classpath.add(configuration.getAmuseHome().getPath() + "/sandbox/ndrost/distributed");

        result.setJavaMain("nl.esciencecenter.amuse.distributed.remote.Pilot");

        List<String> javaArguments = result.getJavaArguments();

        javaArguments.add("--node-label");
        javaArguments.add(nodeLabel);

        javaArguments.add("--resource-name");
        javaArguments.add(resource.getName());

        javaArguments.add("--server-address");
        javaArguments.add(serverAddress);

        String hubs = null;

        if (resource.getHub() != null) {
            hubs = resource.getHub().getAddress();
        } else {
            for (String hub : hubAddresses) {
                if (hubs == null) {
                    hubs = hub;
                } else {
                    hubs = hubs + "," + hub;
                }
            }
        }

        if (hubs != null) {
            javaArguments.add("--hub-addresses");
            javaArguments.add(hubs);
        }

        return result;
    }

    private final int id;

    private final Job job;

    private final Octopus octopus;

    /**
     * @param resource
     * @param queueName
     * @param nodeCount
     * @param timeMinutes
     * @param nodeLabel
     */
    public Reservation(Resource resource, String queueName, int nodeCount, int timeMinutes, String nodeLabel,
            String serverAddress, String[] hubAddresses, Octopus octopus) throws DistributedAmuseException {
        this.octopus = octopus;

        this.id = getNextID();

        try {
            Scheduler scheduler;

            if (resource.isLocal()) {
                scheduler = octopus.jobs().getLocalScheduler();
            } else {
                Credential credential = octopus.credentials().getDefaultCredential("ssh");

                URI uri =
                        new URI(resource.getSchedulerType(), resource.getUsername(), resource.getHostname(), resource.getPort(),
                                null, null, null);
                scheduler = octopus.jobs().newScheduler(uri, credential, null);

            }

            JobDescription jobDescription =
                    createJobDesciption(resource, queueName, nodeCount, timeMinutes, nodeLabel, serverAddress, hubAddresses);
            
            logger.debug("starting job using scheduler " + scheduler);

            this.job = octopus.jobs().submitJob(scheduler, jobDescription);

            logger.debug("submitted job {} with arguments {}", job ,jobDescription.getArguments());

        } catch (Exception e) {
            throw new DistributedAmuseException("cannot start reservation on " + resource.getName(), e);
        }
    }

    public int getID() {
        return id;
    }

    public void cancel() throws DistributedAmuseException {
        try {
            octopus.jobs().cancelJob(job);
        } catch (OctopusIOException | OctopusException e) {
            throw new DistributedAmuseException("failed to cancel job " + job, e);
        }
    }

    /**
     * 
     */
    public void waitUntilStarted() throws DistributedAmuseException {
        try {
            octopus.jobs().waitUntilRunning(job, 0);
        } catch (OctopusIOException | OctopusException e) {
            throw new DistributedAmuseException("failed to get job status " + job, e);
        }
    }

}
