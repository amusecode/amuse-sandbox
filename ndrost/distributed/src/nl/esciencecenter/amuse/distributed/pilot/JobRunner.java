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
package nl.esciencecenter.amuse.distributed.pilot;

import ibis.ipl.Ibis;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.WorkerDescription;
import nl.esciencecenter.amuse.distributed.jobs.JobManager;
import nl.esciencecenter.amuse.distributed.workers.WorkerProxy;

/**
 * A job running on a pilot node.S
 * 
 * @author Niels Drost
 * 
 */
public class JobRunner extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(JobRunner.class);

    private final int jobID;
    private final WorkerProxy workerProxy;
    private final Ibis ibis;

    /**
     * @param jobID
     * @param description
     * @param nodes
     * @param configuration
     * @throws Exception
     */
    public JobRunner(int jobID, WorkerDescription description, AmuseConfiguration configuration, IbisIdentifier[] nodes, Ibis ibis)
            throws Exception {
        this.jobID = jobID;
        this.ibis = ibis;

        File workingDirectory =
                new File(System.getProperty("java.io.tmpdir") + File.separator + "distributed-amuse-worker" + jobID);
        workingDirectory.mkdirs();

        workerProxy = new WorkerProxy(description, configuration, nodes, ibis, workingDirectory);

        setName("Job Runner for " + jobID);
    }

    public void run() {
        try {
            workerProxy.join();
        } catch (InterruptedException e) {
            workerProxy.end();
        }

        //send result message to main node (handled by JobManager, passed to Job)
        try {
            SendPort sendPort = ibis.createSendPort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE);

            IbisIdentifier mainNode = ibis.registry().getElectionResult("amuse");

            sendPort.connect(mainNode, JobManager.PORT_NAME);

            WriteMessage message = sendPort.newMessage();
            
            message.writeInt(jobID);
            
            message.writeString(workerProxy.getResult());

            message.finish();

            sendPort.close();

        } catch (IOException e) {
            logger.error("Failed to report status to main node", e);
        }

    }

}
