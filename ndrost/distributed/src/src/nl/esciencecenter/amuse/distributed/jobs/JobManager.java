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
package nl.esciencecenter.amuse.distributed.jobs;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCreationFailedException;
import ibis.ipl.IbisFactory;
import ibis.ipl.MessageUpcall;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Properties;

import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.WorkerDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class JobManager extends Thread implements MessageUpcall {

    private static final Logger logger = LoggerFactory.getLogger(JobManager.class);

    public static final String PORT_NAME = "job.manager";

    private final Ibis ibis;

    private final ReceivePort receivePort;

    private final PilotNodes nodes;

    //all pending, running, completed, and failed jobs.
    LinkedList<Job> jobs;

    public JobManager(String serverAddress, File tmpDir) throws DistributedAmuseException {
        nodes = new PilotNodes(this);

        try {
            Properties properties = new Properties();
            properties.put("ibis.server.address", serverAddress);
            properties.put("ibis.pool.name", "amuse");
            properties.put("ibis.location", "daemon@local");
            //properties.put("ibis.managementclient", "true");
            //properties.put("ibis.bytescount", "true");

            ibis =
                    IbisFactory.createIbis(DistributedAmuse.IPL_CAPABILITIES, properties, true, nodes,
                            DistributedAmuse.ONE_TO_ONE_PORT_TYPE, DistributedAmuse.MANY_TO_ONE_PORT_TYPE);

            //label this ibis as the master node by running an election with us as the only 
            ibis.registry().elect("amuse");

            ibis.registry().enableEvents();

            receivePort = ibis.createReceivePort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE, PORT_NAME, this);
            receivePort.enableConnections();
            receivePort.enableMessageUpcalls();

        } catch (IOException | IbisCreationFailedException e) {
            throw new DistributedAmuseException("failed to create ibis", e);
        }

        jobs = new LinkedList<Job>();

        //start a thread to run the scheduling
        setName("Job Manager");
        setDaemon(true);
        start();
    }

    public Ibis getIbis() {
        return ibis;
    }

    private synchronized void addJob(Job job) {
        jobs.add(job);

        //run scheduler thread now
        notifyAll();
    }

    /**
     * @return
     */
    private synchronized Job[] getJobs() {
        return jobs.toArray(new Job[jobs.size()]);
    }

    public Job submitWorkerJob(WorkerDescription workerDescription) throws DistributedAmuseException {
        Job result = new Job(workerDescription, ibis);

        addJob(result);

        return result;
    }

    public int submitScriptJob(String script, String arguments, String codeDir, String nodeLabel, boolean useCodeCache)
            throws DistributedAmuseException {
        // TODO Auto-generated method stub
        return 0;
    }


    public int submitPickledJob(String function, String arguments, String nodeLabel) throws DistributedAmuseException {
        // TODO Auto-generated method stub
        return 0;
    }

    private synchronized boolean allBatchJobsDone() {
        for (Job job : jobs) {
            if (job.isBatchJob() && !job.isDone()) {
                return false;
            }
        }
        return true;
    }

    public synchronized void waitForAllBatchJobs() throws DistributedAmuseException {
        while (!allBatchJobsDone()) {
            try {
                wait();
            } catch (InterruptedException e) {
                throw new DistributedAmuseException("Interrupted while waiting for all jobs to finish");
            }
        }
    }
    

    private synchronized Job getJob(int jobID) {
        for (Job job : jobs) {
            if (job.getJobID() == jobID) {
                //will block until the job is finished, then send the result or throw an exception if the job failed, 
                return job;
            }
        }
        return null;
    }

    public synchronized String getJobResult(int jobID) throws DistributedAmuseException {
        for (Job job : jobs) {
            if (job.getJobID() == jobID) {
                //will block until the job is finished, then send the result or throw an exception if the job failed, 
                return job.getResult();
            }
        }
        throw new DistributedAmuseException("Unknown job " + jobID);
    }

    //currently only called by worker connection
    public synchronized void cancelJob(int jobID) throws DistributedAmuseException {
        for (Job job : jobs) {
            if (job.getJobID() == jobID) {
                //blocks while sending a "cancel" message
                job.cancel();
            }
        }
        throw new DistributedAmuseException("Unknown job " + jobID);
    }

    

    public void end() {
        this.interrupt();

        for (Job job : getJobs()) {
            try {
                job.cancel();
            } catch (DistributedAmuseException e) {
                logger.error("Failed to cancel job: " + job, e);
            }
        }
        try {
            logger.debug("Terminating ibis pool");
            ibis.registry().terminate();
        } catch (IOException e) {
            logger.error("Failed to terminate ibis pool", e);
        }
        try {
            ibis.end();
        } catch (IOException e) {
            logger.error("Failed to end ibis", e);
        }
    }

    /**
     * Handles incoming job messages from Pilots
     */
    @Override
    public void upcall(ReadMessage message) throws IOException, ClassNotFoundException {
        int jobID = message.readInt();
        
        logger.debug("Got message for job {}", jobID);
        
        Job job = getJob(jobID);
        
        if (job == null) {
            logger.error("Error handling result for unknown job" + jobID);
            throw new IOException("Error handling result for unknown job" + jobID);
        }
        
        job.handleResult(message);
    }

    /**
     * Wake up the scheduler thread. 
     */
    public synchronized void nudge() {
        notifyAll();
    }

    /**
     * Scheduler thread. Periodically checks if suitable nodes can be found for jobs.
     */
    @Override
    public synchronized void run() {
        while (true) {
            //find nodes for jobs to run on
            Iterator<Job> iterator = jobs.iterator();
            while (iterator.hasNext()) {
                Job job = iterator.next();

                if (job.isPending()) {
                    //find nodes to run this job on. Usually only a single node, but worker jobs may require multiple nodes.
                    PilotNode[] target = nodes.getSuitableNodes(job);

                    //If suitable nodes are found
                    if (target != null) {
                        job.start(target);
                    }
                } else if (job.isObsolete()) {
                    //remove this job
                    iterator.remove();
                }
            }

            try {
                wait(5000);
            } catch (InterruptedException e) {
                logger.debug("Scheduler thread interrupted, time to quit");
                return;
            }
        }
    }

}
