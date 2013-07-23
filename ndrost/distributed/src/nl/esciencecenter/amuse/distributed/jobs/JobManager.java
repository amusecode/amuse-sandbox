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

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.Network;
import nl.esciencecenter.amuse.distributed.WorkerDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class JobManager extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(JobManager.class);

    private final Ibis ibis;

    private final PilotNodes nodes;

    //all pending, running, completed, and failed jobs.
    LinkedList<Job> jobs;

    /**
     * @param tmpDir 
     * @param string
     * @param distributedAmuse
     * @throws DistributedAmuseException
     */
    public JobManager(String serverAddress, File tmpDir) throws DistributedAmuseException {
        nodes = new PilotNodes();

        ibis = Network.createIbis(serverAddress, nodes);
        ibis.registry().enableEvents();

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
     * @param workerDescription
     * @return
     */
    public Job submitWorkerJob(WorkerDescription workerDescription) throws DistributedAmuseException {
        Job result = new Job(workerDescription, ibis);

        addJob(result);

        return result;
    }

    /**
     * @param script
     * @param arguments
     * @param codeDir
     * @param nodeLabel
     * @param useCodeCache
     * @return
     */
    public int submitScriptJob(String script, String arguments, String codeDir, String nodeLabel, boolean useCodeCache)
            throws DistributedAmuseException {
        // TODO Auto-generated method stub
        return 0;
    }

    /**
     * @param function
     * @param arguments
     * @param nodeLabel
     * @return
     */
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

    /**
     * @param jobID
     * @return
     */
    public String getJobResult(int jobID) throws DistributedAmuseException {
        for (Job job : jobs) {
            if (job.getJobID() == jobID) {
                //will block until the job is finished, then send the result or throw an exception if the job failed, 
                return job.getResult();
            }
        }
        throw new DistributedAmuseException("Unknown job " + jobID);
    }

    //currently only called by worker connection
    public void cancelJob(int jobID) throws DistributedAmuseException {
        for (Job job : jobs) {
            if (job.getJobID() == jobID) {
                //blocks while sending a "cancel" message
                job.cancel();
            }
        }
        throw new DistributedAmuseException("Unknown job " + jobID);
    }

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

    public void end() {
        for (Job job : jobs) {
            try {
                job.cancel();
            } catch (DistributedAmuseException e) {
                logger.error("Failed to cancel job: " + job, e);
            }
        }
        try {
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

}
