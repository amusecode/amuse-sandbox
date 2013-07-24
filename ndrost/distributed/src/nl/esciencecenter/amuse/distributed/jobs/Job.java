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
import ibis.ipl.IbisIdentifier;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.util.Arrays;

import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.WorkerDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A job run by the Distributed Amuse system. Contains description and status info, communicates with nodes that actually run job.
 * 
 * @author Niels Drost
 * 
 */
public class Job extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(Job.class);

    private enum State {
        PENDING, INITIALIZING, RUNNING, DONE, REPORTED;
    }

    private static int nextID = 0;

    static synchronized int getNextID() {
        return nextID++;
    }

    private final WorkerDescription workerDescription;

    private final Ibis ibis;

    private final int jobID;

    private State state;

    private PilotNode[] target = null;

    //only for pickled jobs
    private String result = null;

    private Exception error = null;

    //will never timeout until timeout set
    private long timeout = Long.MAX_VALUE;

    public Job(WorkerDescription description, Ibis ibis) {
        this.workerDescription = description;
        this.ibis = ibis;
        this.jobID = getNextID();

        this.state = State.PENDING;
    }

    public WorkerDescription getWorkerDescription() {
        return workerDescription;
    }
    
    public int getNumberOfNodes() {
        if (!isWorkerJob()) {
            return 1;
        }
        return workerDescription.getNrOfNodes();
    }
    
    public String getLabel() {
        //FIXME: support non-worker jobs
        return workerDescription.getNodeLabel();
    }

    public boolean isWorkerJob() {
        return workerDescription != null;
    }

    public boolean isPickledJob() {
        return workerDescription != null;
    }

    public boolean isScriptJob() {
        return workerDescription == null;
    }

    public boolean isBatchJob() {
        return isPickledJob() || isScriptJob();
    }

    public int getJobID() {
        return jobID;
    }

    private synchronized void setState(State newState) {
        state = newState;
        notifyAll();
    }

    public synchronized boolean isPending() {
        return state == State.PENDING;
    }

    public synchronized boolean isRunning() {
        return state == State.RUNNING;
    }

    public synchronized boolean isDone() {
        return state == State.DONE || state == State.REPORTED;
    }

    public synchronized boolean isReported() {
        return state == State.REPORTED;
    }

    public synchronized boolean hasError() {
        return error != null;
    }

    /**
     * Job completely done.
     */
    public synchronized boolean isObsolete() {
        return isReported() && System.currentTimeMillis() > timeout;
    }

    public synchronized void waitUntilRunning() {
        while (!isRunning()) {
            try {
                wait();
            } catch (InterruptedException e) {
                return;
            }
        }
    }

    public synchronized void waitUntilDone() {
        while (!isDone()) {
            try {
                wait();
            } catch (InterruptedException e) {
                return;
            }
        }
    }

    /**
     * @return the result
     * 
     * @throws DistributedAmuseException
     *             in case the job failed for some reason
     */
    public synchronized String getResult() throws DistributedAmuseException {
        if (!isPickledJob()) {
            throw new DistributedAmuseException("Can only get job result for pickled jobs");
        }

        waitUntilDone();

        if (state == State.DONE) {
            setState(State.REPORTED);
        }

        if (hasError()) {
            throw new DistributedAmuseException("Error while running job " + this, error);
        }

        return result;
    }

    /**
     * @param target
     *            the nodes to run this job on.
     */
    public synchronized void start(PilotNode[] target) {
        logger.debug("Running job on target nodes {}", new Object[] {target});
        if (!isPending()) {
            logger.error("Tried to run job {} that was not pending. Ignoring", this);
            return;
        }
        //set state to initializing
        setState(State.INITIALIZING);

        this.target = target;

        //report that this job will be run on the target nodes
        for (PilotNode node : target) {
            node.addJob(this);
        }

        //send out messages to the nodes in a separate thread (see run function below)
        setName("Job " + jobID + " starting thread");
        setDaemon(true);
        start();
    }

    /**
     * @return
     */
    private synchronized IbisIdentifier[] getIbisIdentifiers() {
        if (target == null) {
            return null;
        }

        IbisIdentifier[] result = new IbisIdentifier[target.length];

        for (int i = 0; i < result.length; i++) {
            result[i] = target[i].getIbisIdentifier();
        }

        return result;
    }

    /**
     * Function that starts the job. Only communicates with first node used.
     */
    @Override
    public void run() {
        try {
            PilotNode master = target[0];

            SendPort sendPort = ibis.createSendPort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE);
            ReceivePort receivePort = ibis.createReceivePort(DistributedAmuse.ONE_TO_ONE_PORT_TYPE, null);
            receivePort.enableConnections();

            logger.debug("connecting to pilot");
            
            sendPort.connect(master.getIbisIdentifier(), "pilot");
            
            logger.debug("sending message to pilot");

            WriteMessage writeMessage = sendPort.newMessage();
            
            logger.debug("writing content");

            //where to send the reply
            writeMessage.writeObject(receivePort.identifier());

            //command
            writeMessage.writeString("start");

            //details of job
            writeMessage.writeInt(jobID);
            writeMessage.writeObject(workerDescription);

            writeMessage.writeObject(getIbisIdentifiers());

            //FIXME: transfer files etc

            writeMessage.finish();

            logger.debug("closing sendport");
            
            sendPort.close();
            
            logger.debug("receiving reply from pilot");

            //FIXME: we should use some kind of rpc mechanism
            ReadMessage readMessage = receivePort.receive(60000);

            String statusMessage = readMessage.readString();

            readMessage.finish();
            receivePort.close();

            if (!statusMessage.equals("ok")) {
                setState(State.DONE);
                error = new DistributedAmuseException("Remote node reported error: " + statusMessage);
            }

            setState(State.RUNNING);
            logger.debug("Job {} started on node {}", this, master);
        } catch (IOException e) {
            logger.error("Job failed!", e);
            setState(State.DONE);
            error = e;
        }
    }

    /**
     * 
     */
    public void cancel() throws DistributedAmuseException {
        PilotNode master = target[0];

        try {
            SendPort sendPort = ibis.createSendPort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE);
            ReceivePort receivePort = ibis.createReceivePort(DistributedAmuse.ONE_TO_ONE_PORT_TYPE, null);
            receivePort.enableConnections();

            sendPort.connect(master.getIbisIdentifier(), "jobs");

            WriteMessage writeMessage = sendPort.newMessage();

            //where to send the reply
            writeMessage.writeObject(receivePort.identifier());

            //command
            writeMessage.writeString("cancel");

            writeMessage.writeInt(jobID);

            writeMessage.finish();

            sendPort.close();

            //FIXME: we should use some kind of rpc mechanism
            ReadMessage readMessage = receivePort.receive(60000);

            String statusMessage = readMessage.readString();

            readMessage.finish();
            receivePort.close();

            if (!statusMessage.equals("ok")) {
                setState(State.DONE);
                error = new DistributedAmuseException("Remote node reported error: " + statusMessage);
            }

            setState(State.DONE);
            logger.debug("Job {} canceled on node {}", this, master);

        } catch (IOException e) {
            throw new DistributedAmuseException("Failed to cancel job " + this, e);
        }
    }

    void handleResultMessage(ReadMessage message) throws DistributedAmuseException {
        try {
            String statusMessage = message.readString();

            if (!statusMessage.equals("ok")) {
                logger.warn("Job ended in error: " + statusMessage);
            }

            this.error = (Exception) message.readObject();

            //FIXME: read result files

            setState(State.DONE);
            logger.debug("Job {} done on node {}", this, message.origin());
        } catch (IOException | ClassNotFoundException e) {
            throw new DistributedAmuseException("Failed to handle job result " + this, e);
        }
    }

    @Override
    public String toString() {
        return "Job [jobID=" + jobID + ", label=" + getLabel() + ", state=" + state + ", target=" + Arrays.toString(target) + ", result=" + result
                + ", error=" + error + ", timeout=" + timeout + "]";
    }

    


}
